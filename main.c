#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>

#define MAX_SYMBOLS 256 // Number of possible byte values (0-255)

// Node structure for Huffman Tree
typedef struct Node {
    unsigned char symbol;          // Symbol represented by this node
    unsigned int frequency;        // Frequency of this symbol
    struct Node* left;             // Pointer to the left child
    struct Node* right;            // Pointer to the right child
} Node;

// Priority Queue for building the Huffman Tree
typedef struct PriorityQueue {
    Node* nodes[MAX_SYMBOLS];
    int size;
} PriorityQueue;

// Table for Huffman encoding
typedef struct {
    char codes[MAX_SYMBOLS][MAX_SYMBOLS]; // Max code length = 256 bits
    int lengths[MAX_SYMBOLS];             // Lengths of the codes for each symbol
} HuffmanTable;

// Function prototypes
Node* build_huffman_tree(unsigned int frequencies[]);
void generate_huffman_table(Node* root, HuffmanTable* table, char* code, int length);
void compress_file(const char* input_file, const char* output_file);
void decompress_file(const char* input_file, const char* output_file);
void free_huffman_tree(Node* root);
void show_progress(unsigned long long processed, unsigned long long total, const char* action);

PriorityQueue* create_priority_queue();
void priority_queue_push(PriorityQueue* pq, Node* node);
Node* priority_queue_pop(PriorityQueue* pq);

// I/O Bit-level functions
void write_bit(FILE* file, int bit, unsigned char* buffer, int* bit_pos);
void flush_bits(FILE* file, unsigned char* buffer, int* bit_pos);
int read_bit(FILE* file, unsigned char* buffer, int* bit_pos, int* buffer_size);

// Priority queue functions
PriorityQueue* create_priority_queue() {
    PriorityQueue* pq = (PriorityQueue*)malloc(sizeof(PriorityQueue));
    pq->size = 0;
    return pq;
}

void priority_queue_push(PriorityQueue* pq, Node* node) {
    pq->nodes[pq->size++] = node;

    for (int i = pq->size - 1; i > 0; i--) {
        int parent = (i - 1) / 2;
        if (pq->nodes[parent]->frequency > pq->nodes[i]->frequency) {
            Node* tmp = pq->nodes[parent];
            pq->nodes[parent] = pq->nodes[i];
            pq->nodes[i] = tmp;
        } else {
            break;
        }
    }
}

Node* priority_queue_pop(PriorityQueue* pq) {
    Node* top = pq->nodes[0];
    pq->nodes[0] = pq->nodes[--pq->size];

    int i = 0;
    while (1) {
        int left = 2 * i + 1;
        int right = 2 * i + 2;
        int smallest = i;

        if (left < pq->size && pq->nodes[left]->frequency < pq->nodes[smallest]->frequency)
            smallest = left;
        if (right < pq->size && pq->nodes[right]->frequency < pq->nodes[smallest]->frequency)
            smallest = right;

        if (smallest != i) {
            Node* tmp = pq->nodes[i];
            pq->nodes[i] = pq->nodes[smallest];
            pq->nodes[smallest] = tmp;
            i = smallest;
        } else {
            break;
        }
    }
    return top;
}

// Build Huffman tree based on frequencies
Node* build_huffman_tree(unsigned int frequencies[]) {
    PriorityQueue* pq = create_priority_queue();

    for (int i = 0; i < MAX_SYMBOLS; i++) {
        if (frequencies[i] > 0) {
            Node* node = (Node*)malloc(sizeof(Node));
            node->symbol = i;
            node->frequency = frequencies[i];
            node->left = node->right = NULL;
            priority_queue_push(pq, node);
        }
    }

    while (pq->size > 1) {
        Node* left = priority_queue_pop(pq);
        Node* right = priority_queue_pop(pq);

        Node* parent = (Node*)malloc(sizeof(Node));
        parent->symbol = 0;  // Internal nodes do not represent symbols
        parent->frequency = left->frequency + right->frequency;
        parent->left = left;
        parent->right = right;

        priority_queue_push(pq, parent);
    }

    Node* root = priority_queue_pop(pq);
    free(pq);
    return root;
}

// Generate Huffman encoding table
void generate_huffman_table(Node* root, HuffmanTable* table, char* code, int length) {
    if (!root->left && !root->right) {
        code[length] = '\0';
        strcpy(table->codes[root->symbol], code);
        table->lengths[root->symbol] = length;
        return;
    }

    if (root->left) {
        code[length] = '0';
        generate_huffman_table(root->left, table, code, length + 1);
    }

    if (root->right) {
        code[length] = '1';
        generate_huffman_table(root->right, table, code, length + 1);
    }
}

// Display progress
void show_progress(unsigned long long processed, unsigned long long total, const char* action) {
    int progress = (int)((processed * 100) / total);
    printf("\r%s: %d%%", action, progress);
    fflush(stdout);
}

// Write a bit to the output file
void write_bit(FILE* file, int bit, unsigned char* buffer, int* bit_pos) {
    if (bit)
        *buffer |= (1 << *bit_pos);

    if (++(*bit_pos) == 8) {
        fwrite(buffer, 1, 1, file);
        *buffer = 0;
        *bit_pos = 0;
    }
}

// Flush remaining bits
void flush_bits(FILE* file, unsigned char* buffer, int* bit_pos) {
    if (*bit_pos > 0) {
        fwrite(buffer, 1, 1, file);
        *bit_pos = 0;
    }
}

// Read a bit from the input file
int read_bit(FILE* file, unsigned char* buffer, int* bit_pos, int* buffer_size) {
    if (*bit_pos == 0) {
        *buffer_size = fread(buffer, 1, 1, file);
        if (*buffer_size == 0)
            return -1;
    }

    int bit = (*buffer >> *bit_pos) & 1;
    *bit_pos = (*bit_pos + 1) % 8;
    return bit;
}

// Compress a file
void compress_file(const char* input_file, const char* output_file) {
    FILE* in = fopen(input_file, "rb");
    FILE* out = fopen(output_file, "wb");

    if (!in || !out) {
        fprintf(stderr, "Error: Could not open files.\n");
        return;
    }

    // Build frequency table
    unsigned int frequencies[MAX_SYMBOLS] = {0};
    unsigned char buffer[1];
    unsigned long long input_size = 0;

    while (fread(buffer, 1, 1, in)) {
        frequencies[buffer[0]]++;
        input_size++;
    }

    // Write input file size to the output file
    fwrite(&input_size, sizeof(input_size), 1, out);

    // Build Huffman tree
    Node* root = build_huffman_tree(frequencies);
    HuffmanTable table = {0};
    char code[MAX_SYMBOLS];
    generate_huffman_table(root, &table, code, 0);

    // Write frequency table to the output file
    fwrite(frequencies, sizeof(frequencies), 1, out);

    // Write encoded data
    rewind(in);
    unsigned char bit_buffer = 0;
    int bit_pos = 0;
    unsigned long long bytes_processed = 0;

    while (fread(buffer, 1, 1, in)) {
        char* code = table.codes[buffer[0]];
        int length = table.lengths[buffer[0]];

        for (int i = 0; i < length; i++) {
            write_bit(out, code[i] - '0', &bit_buffer, &bit_pos);
        }

        bytes_processed++;
        if (bytes_processed % 1024 == 0 || bytes_processed == input_size) {
            show_progress(bytes_processed, input_size, "Compressing");
        }
    }

    flush_bits(out, &bit_buffer, &bit_pos);
    printf("\n");
    free_huffman_tree(root);

    fclose(in);
    fclose(out);
}

// Decompress a file
void decompress_file(const char* input_file, const char* output_file) {
    FILE* in = fopen(input_file, "rb");
    FILE* out = fopen(output_file, "wb");

    if (!in || !out) {
        fprintf(stderr, "Error: Could not open files.\n");
        return;
    }

    // Read input file size
    unsigned long long input_size;
    fread(&input_size, sizeof(input_size), 1, in);

    // Read frequency table
    unsigned int frequencies[MAX_SYMBOLS];
    fread(frequencies, sizeof(frequencies), 1, in);

    // Rebuild Huffman tree
    Node* root = build_huffman_tree(frequencies);

    // Decode data
    Node* current = root;
    unsigned char bit_buffer = 0;
    int bit_pos = 0;
    int buffer_size = 0;

    unsigned long long decoded = 0;
    while (decoded < input_size) {
        int bit = read_bit(in, &bit_buffer, &bit_pos, &buffer_size);
        if (bit == -1) // EOF
            break;

        current = bit ? current->right : current->left;

        if (!current->left && !current->right) {
            fwrite(&current->symbol, 1, 1, out);
            current = root;
            decoded++;

            if (decoded % 1024 == 0 || decoded == input_size) {
                show_progress(decoded, input_size, "Decompressing");
            }
        }
    }

    printf("\n");
    free_huffman_tree(root);
    fclose(in);
    fclose(out);
}

// Free the Huffman tree
void free_huffman_tree(Node* root) {
    if (root) {
        free_huffman_tree(root->left);
        free_huffman_tree(root->right);
        free(root);
    }
}

// Main function
int main(int argc, char** argv) {
    if (argc != 4) {
        printf("Usage: %s compress|decompress input_file output_file\n", argv[0]);
        return 1;
    }

    if (strcmp(argv[1], "compress") == 0) {
        compress_file(argv[2], argv[3]);
    } else if (strcmp(argv[1], "decompress") == 0) {
        decompress_file(argv[2], argv[3]);
    } else {
        printf("Unknown command: %s\n", argv[1]);
    }

    return 0;
}