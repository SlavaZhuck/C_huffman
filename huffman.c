#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>

/* Structure for Huffman tree node */
typedef struct HuffmanNode {
    unsigned char data;           // Character
    unsigned int frequency;       // Frequency of the character
    struct HuffmanNode *left;     // Left child
    struct HuffmanNode *right;    // Right child
} HuffmanNode;

/* Structure for Min Heap Node */
typedef struct MinHeapNode {
    HuffmanNode **array;          // Array of huffman nodes
    unsigned int size;            // Current size of min heap
    unsigned int capacity;        // Capacity of min heap
} MinHeap;

/* Function to allocate a new Huffman tree node */
HuffmanNode* newHuffmanNode(unsigned char data, unsigned int freq);

/* Function to create a min heap of given capacity */
MinHeap* createMinHeap(unsigned int capacity);

/* Function to swap two min heap nodes */
void swapMinHeapNodes(HuffmanNode** a, HuffmanNode** b);

/* Standard function to heapify at a given index */
void minHeapify(MinHeap* minHeap, int index);

/* Function to check if size of heap is 1 */
bool isSizeOne(MinHeap* minHeap);

/* Function to extract minimum value node from heap */
HuffmanNode* extractMin(MinHeap* minHeap);

/* Function to insert a new node to Min Heap */
void insertMinHeap(MinHeap* minHeap, HuffmanNode* minHeapNode);

/* Function to build the min heap */
void buildMinHeap(MinHeap* minHeap);

/* Function to print huffman codes from the root of Huffman Tree */
void printCodes(HuffmanNode* root, int arr[], int top);

/* Function to build Huffman Tree and print codes */
HuffmanNode* buildHuffmanTree(unsigned char data[], int freq[], int size);

/* Main function for Huffman codes using min heap */
void HuffmanCodes(unsigned char data[], int freq[], int size);

/* Function to compress a file using Huffman coding */
void compressFile(const char* input_file, const char* output_file);

/* Function to decompress a file compressed with Huffman coding */
void decompressFile(const char* input_file, const char* output_file);

/* Main implementation functions (to be implemented) */
// void buildFrequencyTable(const char* filename, unsigned char* data, int* freq, int* size);
// void writeCompressedFile(const char* output_file, HuffmanNode* root, const char* input_file);
// HuffmanNode* readHuffmanTree(FILE* input);
// void writeHuffmanTree(FILE* output, HuffmanNode* root);

/* Function to allocate a new Huffman tree node */
HuffmanNode* newHuffmanNode(unsigned char data, unsigned int freq) {
    HuffmanNode* node = (HuffmanNode*)malloc(sizeof(HuffmanNode));
    if (node == NULL) {
        fprintf(stderr, "Memory allocation failed\n");
        exit(1);
    }
    
    node->data = data;
    node->frequency = freq;
    node->left = NULL;
    node->right = NULL;
    
    return node;
}

/* Function to create a min heap of given capacity */
MinHeap* createMinHeap(unsigned int capacity) {
    MinHeap* minHeap = (MinHeap*)malloc(sizeof(MinHeap));
    if (minHeap == NULL) {
        fprintf(stderr, "Memory allocation failed\n");
        exit(1);
    }
    
    // Initialize size as 0
    minHeap->size = 0;
    minHeap->capacity = capacity;
    
    // Allocate memory for array of nodes
    minHeap->array = (HuffmanNode**)malloc(capacity * sizeof(HuffmanNode*));
    if (minHeap->array == NULL) {
        fprintf(stderr, "Memory allocation failed\n");
        free(minHeap);
        exit(1);
    }
    
    return minHeap;
}

/* Function to swap two min heap nodes */
void swapMinHeapNodes(HuffmanNode** a, HuffmanNode** b) {
    HuffmanNode* temp = *a;
    *a = *b;
    *b = temp;
}

/* Standard function to heapify at a given index */
void minHeapify(MinHeap* minHeap, int idx) {
    int smallest = idx;
    int left = 2 * idx + 1;
    int right = 2 * idx + 2;
    
    // If left child is smaller than root
    if (left < minHeap->size && 
        minHeap->array[left]->frequency < minHeap->array[smallest]->frequency) {
        smallest = left;
    }
    
    // If right child is smaller than smallest so far
    if (right < minHeap->size && 
        minHeap->array[right]->frequency < minHeap->array[smallest]->frequency) {
        smallest = right;
    }
    
    // If smallest is not root
    if (smallest != idx) {
        swapMinHeapNodes(&minHeap->array[smallest], &minHeap->array[idx]);
        minHeapify(minHeap, smallest);
    }
}

/* Function to check if size of heap is 1 */
bool isSizeOne(MinHeap* minHeap) {
    return (minHeap->size == 1);
}

/* Function to extract minimum value node from heap */
HuffmanNode* extractMin(MinHeap* minHeap) {
    if (minHeap->size == 0) {
        return NULL;
    }
    
    HuffmanNode* temp = minHeap->array[0];
    
    // Replace the root with the last element
    minHeap->array[0] = minHeap->array[minHeap->size - 1];
    
    // Reduce heap size
    minHeap->size--;
    
    // Heapify the root
    minHeapify(minHeap, 0);
    
    return temp;
}

/* Function to insert a new node to Min Heap */
void insertMinHeap(MinHeap* minHeap, HuffmanNode* minHeapNode) {
    if (minHeap->size >= minHeap->capacity) {
        fprintf(stderr, "Min Heap overflow\n");
        return;
    }
    
    // First insert the new node at the end
    minHeap->size++;
    int i = minHeap->size - 1;
    minHeap->array[i] = minHeapNode;
    
    // Fix the min heap property if it is violated
    // i.e., child's frequency is less than parent's frequency
    while (i && minHeap->array[i]->frequency < 
                 minHeap->array[(i - 1) / 2]->frequency) {
        swapMinHeapNodes(&minHeap->array[i], &minHeap->array[(i - 1) / 2]);
        i = (i - 1) / 2;
    }
}

/* Function to build the min heap */
void buildMinHeap(MinHeap* minHeap) {
    int n = minHeap->size - 1;
    for (int i = (n - 1) / 2; i >= 0; --i) {
        minHeapify(minHeap, i);
    }
}

/* Function to build Huffman tree */
HuffmanNode* buildHuffmanTree(unsigned char data[], int freq[], int size) {
    // Create a min heap of capacity equal to size
    MinHeap* minHeap = createMinHeap(size);
    
    // For each character, create a huffman node and insert it into min heap
    for (int i = 0; i < size; ++i) {
        minHeap->array[i] = newHuffmanNode(data[i], freq[i]);
    }
    
    // Set size of min heap
    minHeap->size = size;
    
    // Build the min heap
    buildMinHeap(minHeap);
    
    // Iterate while size of heap doesn't become 1
    while (!isSizeOne(minHeap)) {
        // Extract the two minimum frequency nodes
        HuffmanNode* left = extractMin(minHeap);
        HuffmanNode* right = extractMin(minHeap);
        
        // Create a new internal node with frequency equal to the sum of the
        // two extracted nodes. Make the extracted nodes as left and right children
        // of this new node. '$' is a special character used for internal nodes
        HuffmanNode* top = newHuffmanNode('$', left->frequency + right->frequency);
        top->left = left;
        top->right = right;
        
        // Add this node to the min heap
        insertMinHeap(minHeap, top);
    }
    
    // The remaining node is the root node and the tree is complete
    HuffmanNode* root = extractMin(minHeap);
    
    // Free the min heap but not the tree nodes
    free(minHeap->array);
    free(minHeap);
    
    return root;
}

/* Function to print huffman codes from the root of Huffman Tree.
 * It uses arr[] to store codes */
void printCodes(HuffmanNode* root, int arr[], int top) {
    // Assign 0 to left edge and 1 to right edge
    if (root->left) {
        arr[top] = 0;
        printCodes(root->left, arr, top + 1);
    }
    
    if (root->right) {
        arr[top] = 1;
        printCodes(root->right, arr, top + 1);
    }
    
    // If this is a leaf node, then it contains a character
    // Print the character and its code
    if (!root->left && !root->right) {
        printf("%c: ", root->data);
        
        // Print the code
        for (int i = 0; i < top; ++i) {
            printf("%d", arr[i]);
        }
        printf("\n");
    }
}

/* Main function that builds a Huffman Tree and print codes */
void HuffmanCodes(unsigned char data[], int freq[], int size) {
    // Construct Huffman Tree
    HuffmanNode* root = buildHuffmanTree(data, freq, size);
    
    // Print Huffman codes using the Huffman tree built above
    int arr[100], top = 0; // Assuming maximum code length is 100
    
    printf("Huffman Codes:\n");
    printCodes(root, arr, top);
}

/* Function to build frequency table from input file */
void buildFrequencyTable(const char* filename, unsigned char* data, int* freq, int* size) {
    FILE* file = fopen(filename, "rb");
    if (file == NULL) {
        fprintf(stderr, "Could not open file %s\n", filename);
        exit(1);
    }
    
    // Initialize frequency array to 0
    for (int i = 0; i < 256; i++) {
        freq[i] = 0;
    }
    
    // Count frequency of each character
    unsigned char ch;
    *size = 0;
    
    while (fread(&ch, 1, 1, file) == 1) {
        if (freq[ch] == 0) {
            data[*size] = ch;
            (*size)++;
        }
        freq[ch]++;
    }
    
    fclose(file);
}

/* Function to write a bit to output file */
void writeBit(FILE* output, int bit, unsigned char* buffer, int* bitCount) {
    // Set the bit in the buffer
    if (bit) {
        *buffer |= (1 << (7 - *bitCount));
    }
    
    (*bitCount)++;
    
    // If buffer is full, write it to the file
    if (*bitCount == 8) {
        fwrite(buffer, 1, 1, output);
        *buffer = 0;
        *bitCount = 0;
    }
}

/* Function to write Huffman tree to output file for decompression */
void writeHuffmanTree(FILE* output, HuffmanNode* root, unsigned char* buffer, int* bitCount) {
    // If internal node, write 0 and recursively write left and right subtrees
    if (root->left || root->right) {
        writeBit(output, 0, buffer, bitCount);
        writeHuffmanTree(output, root->left, buffer, bitCount);
        writeHuffmanTree(output, root->right, buffer, bitCount);
    } 
    // If leaf node, write 1 followed by the character
    else {
        writeBit(output, 1, buffer, bitCount);
        
        // Write the 8 bits of the character
        for (int i = 0; i < 8; i++) {
            int bit = (root->data >> (7 - i)) & 1;
            writeBit(output, bit, buffer, bitCount);
        }
    }
}

/* Function to generate Huffman codes and store them in a table */
void generateCodes(HuffmanNode* root, char* code, int depth, char** codes) {
    if (root == NULL) return;
    
    // If this is a leaf node, store the code
    if (!root->left && !root->right) {
        codes[root->data] = strdup(code);
        return;
    }
    
    // Traverse left
    code[depth] = '0';
    code[depth + 1] = '\0';
    generateCodes(root->left, code, depth + 1, codes);
    
    // Traverse right
    code[depth] = '1';
    code[depth + 1] = '\0';
    generateCodes(root->right, code, depth + 1, codes);
}

/* Function to compress a file using Huffman coding */
void compressFile(const char* input_file, const char* output_file) {
    unsigned char data[256];
    int freq[256];
    int size = 0;
    
    // Build frequency table
    buildFrequencyTable(input_file, data, freq, &size);
    
    // Build Huffman tree
    HuffmanNode* root = buildHuffmanTree(data, freq, size);
    
    // Generate Huffman codes
    char* codes[256] = {NULL};
    char code[256] = {0};
    generateCodes(root, code, 0, codes);
    
    // Open input and output files
    FILE* input = fopen(input_file, "rb");
    FILE* output = fopen(output_file, "wb");
    
    if (input == NULL || output == NULL) {
        fprintf(stderr, "Error opening files\n");
        exit(1);
    }
    
    // Write the size of the frequency table
    fwrite(&size, sizeof(int), 1, output);
    
    // Write the Huffman tree to output file
    unsigned char buffer = 0;
    int bitCount = 0;
    writeHuffmanTree(output, root, &buffer, &bitCount);
    
    // If there are any bits remaining in the buffer, write them
    if (bitCount > 0) {
        fwrite(&buffer, 1, 1, output);
        buffer = 0;
        bitCount = 0;
    }
    
    // Compress data
    unsigned char ch;
    
    // Get file size for progress reporting
    fseek(input, 0, SEEK_END);
    long fileSize = ftell(input);
    fseek(input, 0, SEEK_SET);
    
    printf("Compressing file...\n");
    long bytesRead = 0;
    
    while (fread(&ch, 1, 1, input) == 1) {
        char* code_str = codes[ch];
        if (code_str != NULL) {
            for (int i = 0; code_str[i]; i++) {
                int bit = code_str[i] - '0';
                writeBit(output, bit, &buffer, &bitCount);
            }
        }
        
        bytesRead++;
        if (bytesRead % (fileSize / 10) == 0) {
            printf("Progress: %.0f%%\n", (double)bytesRead / fileSize * 100);
        }
    }
    
    // If there are any bits remaining in the buffer, write them
    if (bitCount > 0) {
        fwrite(&buffer, 1, 1, output);
    }
    
    // Close files
    fclose(input);
    fclose(output);
    
    // Free memory
    for (int i = 0; i < 256; i++) {
        if (codes[i] != NULL) {
            free(codes[i]);
        }
    }
    
    printf("Compression complete. File saved as %s\n", output_file);
}

/* Function to read a bit from input file */
int readBit(FILE* input, unsigned char* buffer, int* bitCount) {
    // If buffer is empty, read a new byte
    if (*bitCount == 0) {
        if (fread(buffer, 1, 1, input) != 1) {
            return -1; // Error or EOF
        }
        *bitCount = 8;
    }
    
    // Extract the next bit from the buffer
    int bit = (*buffer >> (8 - *bitCount)) & 1;
    (*bitCount)--;
    
    return bit;
}

/* Function to read and reconstruct the Huffman tree from a compressed file */
HuffmanNode* readHuffmanTree(FILE* input, unsigned char* buffer, int* bitCount) {
    // Read a bit to determine if it's a leaf node or internal node
    int bit = readBit(input, buffer, bitCount);
    if (bit == -1) return NULL; // Error or EOF
    
    // If it's a leaf node (bit is 1)
    if (bit == 1) {
        // Read the 8 bits of the character
        unsigned char ch = 0;
        for (int i = 0; i < 8; i++) {
            bit = readBit(input, buffer, bitCount);
            if (bit == -1) return NULL; // Error or EOF
            
            ch = (ch << 1) | bit;
        }
        
        // Create and return a new leaf node
        return newHuffmanNode(ch, 0); // Frequency is not needed for decompression
    }
    // If it's an internal node (bit is 0)
    else {
        // Create a new internal node
        HuffmanNode* node = newHuffmanNode('$', 0); // '$' and 0 are placeholder values
        
        // Recursively read left and right children
        node->left = readHuffmanTree(input, buffer, bitCount);
        node->right = readHuffmanTree(input, buffer, bitCount);
        
        return node;
    }
}

/* Function to decompress a file that was compressed using Huffman compression */
void decompressFile(const char* input_file, const char* output_file) {
    FILE* input = fopen(input_file, "rb");
    FILE* output = fopen(output_file, "wb");
    
    if (input == NULL || output == NULL) {
        fprintf(stderr, "Error opening files\n");
        exit(1);
    }
    
    // Read the size of the frequency table
    int size;
    if (fread(&size, sizeof(int), 1, input) != 1) {
        fprintf(stderr, "Error reading file size\n");
        exit(1);
    }
    
    printf("Decompressing file...\n");
    
    // Read and reconstruct the Huffman tree
    unsigned char buffer = 0;
    int bitCount = 0;
    HuffmanNode* root = readHuffmanTree(input, &buffer, &bitCount);
    
    if (root == NULL) {
        fprintf(stderr, "Error reading Huffman tree\n");
        exit(1);
    }
    
    // Decompress data
    HuffmanNode* current = root;
    int bit;
    
    // Keep reading bits until we reach EOF
    while ((bit = readBit(input, &buffer, &bitCount)) != -1) {
        // Navigate the Huffman tree based on the bit
        if (bit == 0) {
            current = current->left;
        } else {
            current = current->right;
        }
        
        // If we've reached a leaf node, write the character and reset to the root
        if (current->left == NULL && current->right == NULL) {
            fwrite(&current->data, 1, 1, output);
            current = root;
        }
    }
    
    // Close files
    fclose(input);
    fclose(output);
    
    // Free the Huffman tree
    freeHuffmanTree(root);
    
    printf("Decompression complete. File saved as %s\n", output_file);
}

/* Function to free Huffman tree recursively */
void freeHuffmanTree(HuffmanNode* node) {
    if (node == NULL) return;
    
    freeHuffmanTree(node->left);
    freeHuffmanTree(node->right);
    free(node);
}

int main(int argc, char* argv[]) {
    if (argc != 4) {
        printf("Usage: %s <mode: c|d> <input_file> <output_file>\n", argv[0]);
        printf("  c: compress\n");
        printf("  d: decompress\n");
        return 1;
    }
    
    char mode = argv[1][0];
    const char* input_file = argv[2];
    const char* output_file = argv[3];
    
    printf("Input file: %s\n", input_file);
    printf("Output file: %s\n", output_file);
    
    if (mode == 'c') {
        printf("Huffman compression starting...\n");
        // Compress file
        compressFile(input_file, output_file);
    } else if (mode == 'd') {
        printf("Huffman decompression starting...\n");
        // Decompress file
        decompressFile(input_file, output_file);
    } else {
        printf("Invalid mode. Use 'c' for compression or 'd' for decompression.\n");
        return 1;
    }
    
    // For demonstration, if we're compressing, print the Huffman codes of a small example
    if (mode == 'c') {
        printf("\nDemo of Huffman coding on a small example:\n");
        unsigned char data[] = {'a', 'b', 'c', 'd', 'e', 'f'};
        int freq[] = {5, 9, 12, 13, 16, 45};
        int size = sizeof(data) / sizeof(data[0]);
        
        HuffmanCodes(data, freq, size);
    }
    
    return 0;
}
