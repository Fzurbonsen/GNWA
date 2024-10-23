/*

    gnwa.c
    This file contains the implementation of the Graph-Needleman-Wunsch-Algorithm.
    Author: Frederic zur Bonsen <fzurbonsen@student.ethz.ch>

*/

#include "gnwa.h"


// Function to find the max of three values
int32_t gnwa_max(int32_t a, int32_t b, int32_t c) {
    if (a >= b && a >= c) return a;
    if (b >= a && b >= c) return b;
    return c;
}


// Function to transform char seqeunce into num sequence
int8_t* gnwa_create_num(const char* seq,
                        const int32_t len,
                        const int8_t* nt_table) {

    int8_t* num = (int8_t*)malloc(len);
    for (int i = 0; i < len; ++i) num[i] = nt_table[(int)seq[i]];
    return num;
}


// Function to print a num sequence
void gnwa_num_print(FILE* file, int8_t* num, int32_t len) {
    for (int i = 0; i < len; ++i) {
        fprintf(file, "%i", num[i]);
    }
    fprintf(file, "\n");
}


// Function to create a score matrix
int8_t* gnwa_create_score_matrix(int32_t match, int32_t mismatch) {
    int8_t* mat = (int8_t*)calloc(25, sizeof(int8_t));
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
            mat[i * 4 + j] = i == j ? match : -mismatch;
        }
    }
    return mat;
}


// Function to create nt table
int8_t* gnwa_create_nt_table(void) {
    int8_t* ret_nt_table = calloc(128, sizeof(int8_t));
    int8_t nt_table[128] = {
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  3, 0, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  3, 0, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4
    };
    memcpy(ret_nt_table, nt_table, 128*sizeof(int8_t));
    return ret_nt_table;
}


// Function to print a CIGAR
void gnwa_cigar_print(FILE* file, gnwa_cigar_t* cigar) {
    for (int i = 0; i < cigar->length; ++i) {
        fprintf(file, "%i%c", cigar->elements[i].length, cigar->elements[i].type);
    }
    fprintf(file, "\n");
}


// Function to destroy a CIGAR
void gnwa_cigar_destroy(gnwa_cigar_t* cigar) {
    free(cigar->elements);
    cigar->elements = NULL;
    cigar->length = 0;
}


// Function to create a new node struct
gnwa_node_t* gnwa_node_create(const uint64_t id,
                                const char* seq,
                                int8_t* nt_table,
                                int8_t* score_matrix) {

    gnwa_node_t* node = calloc(1, sizeof(gnwa_node_t));

    // Prepare memory for the different fields
    node->seq = (char*)malloc(strlen(seq) + 1);
    node->num = (int8_t*)malloc(strlen(seq));

    // Fill the fields with the respective values
    node->id = id;
    node->len = strlen(seq);
    node->count_prev = 0;
    node->count_next = 0;
    strcpy(node->seq, seq); // copy string
    node->num = gnwa_create_num(seq, strlen(seq), nt_table); // transform char to int for sequence

    return node;
}


// Function to destroy a node
void gnwa_node_destroy(gnwa_node_t* node) {
    free(node->seq);
    free(node->num);
    free(node->prev);
    free(node->next);
    free(node);
}


// Funtion to add a previous node
void gnwa_node_add_prev(gnwa_node_t* node, gnwa_node_t* prev) {
    node->count_prev++;
    node->prev = (gnwa_node_t**)realloc(node->prev, node->count_prev * sizeof(gnwa_node_t*));
    node->prev[node->count_prev - 1] = prev;
}


// Function to add a next node
void gnwa_node_add_next(gnwa_node_t* node, gnwa_node_t* next) {
    node->count_next++;
    node->next = (gnwa_node_t**)realloc(node->next, node->count_next * sizeof(gnwa_node_t*));
    node->next[node->count_next - 1] = next;
}


// Function to delete a previous node
void gnwa_node_del_prev(gnwa_node_t* node, gnwa_node_t* prev) {
    gnwa_node_t** new_prev = (gnwa_node_t**)malloc(node->count_prev * sizeof(gnwa_node_t*));
    gnwa_node_t** iter = node->prev;

    for (int i = 0; i < node->count_prev; ++i, ++iter) {
        if (*iter != prev) {
            new_prev[i] = *iter;
        }
    }
    free(node->prev);
    node->prev = new_prev;
    node->count_prev--;
}


// Function to delete a next node
void gnwa_node_del_next(gnwa_node_t* node, gnwa_node_t* next) {
    gnwa_node_t** new_next = (gnwa_node_t**)malloc(node->count_next * sizeof(gnwa_node_t*));
    gnwa_node_t** iter = node->next;

    for (int i = 0; i < node->count_next; ++i, ++iter) {
        if (*iter != next) {
            new_next[i] = *iter;
        }
    }
    free(node->next);
    node->next = new_next;
    node->count_next--;
}


// Function to add an edge between two nodes
void gnwa_node_add_edge(gnwa_node_t* start, gnwa_node_t* end) {
    // Check if the edge already exists
    for (int i = 0; i < start->count_next; ++i) {
        if (start->next[i] == end) return;
    }
    gnwa_node_add_next(start, end);
    gnwa_node_add_prev(end, start);
}


// Function to delete an edge between two nodes
void gnwa_node_del_edge(gnwa_node_t* start, gnwa_node_t* end) {
    gnwa_node_del_next(start, end);
    gnwa_node_del_prev(end, start);
}


// Function to create a path
gnwa_path_t* gnwa_path_create(gnwa_node_t** nodes, uint32_t len) {
    // Allocate memory
    gnwa_path_t* path = (gnwa_path_t*)malloc(sizeof(gnwa_path_t));
    path->nodes = (gnwa_node_t**)malloc(len * sizeof(gnwa_node_t*));
    // Fill the fields
    path->len = len;
    for (int i = 0; i < len; ++i) {
        path->nodes[i] = nodes[i];
    }
    return path;
}


// Function to destroy a path
void gnwa_path_destroy(gnwa_path_t* path) {
    free(path->nodes);
    free(path);
}


// Function to get the complete sequence of a path
char* gnwa_path_get_sequence(gnwa_path_t* path) {
    int32_t len = 0;
    // Get the length of the sequence
    for (int i = 0; i < path->len; ++i) {
        len += path->nodes[i]->len;
    }

    // Allocate memory for the sequence
    char* seq = (char*)malloc(len * sizeof(char) + 1);

    // Iterate over all nodes to build sequence
    int32_t pos = 0;
    for (int i = 0; i < path->len; ++i) {
        strcpy(seq + pos, path->nodes[i]->seq);
        pos += path->nodes[i]->len; 
    }

    seq[len] = '\0';
    return seq;
}


// Function to get the complete num sequence of a path
int32_t gnwa_path_get_num_sequence(gnwa_path_t* path, int8_t** num) {
    int32_t len = 0;
    
    // Get the total length of the sequence from all nodes
    for (int i = 0; i < path->len; ++i) {
        len += path->nodes[i]->len;
    }

    // Reallocate memory only if num is NULL, otherwise reuse existing memory
    if (*num != NULL) {
        free(*num);
    }
    *num = (int8_t*)malloc(len * sizeof(int8_t));

    // Iterate over all nodes and copy their numerical sequences into num
    int32_t pos = 0;
    for (int i = 0; i < path->len; ++i) {
        for (int j = 0; j < path->nodes[i]->len; ++j, ++pos) {
            (*num)[pos] = path->nodes[i]->num[j];
        }
    }

    return len;
}



// Function to print a path
void gnwa_path_print(FILE* file, gnwa_path_t* path) {
    for (int i = 0; i < path->len; ++i) {
        fprintf(file, "[%li]->", path->nodes[i]->id);
    }
    fprintf(file, "\n");
}


// Function to create a graph
gnwa_graph_t* gnwa_graph_create(uint32_t size) {
    gnwa_graph_t* graph = calloc(1, sizeof(gnwa_graph_t));
    graph->nodes = malloc(size * sizeof(gnwa_node_t*));
    if (!graph || !graph->nodes) {
        fprintf(stderr,  "Error:[gnwa] could not allocate memory for the graph of %u nodes!\n", size);
        exit(1);
    }
    return graph;
}


// Function to find the max node of a graph
gnwa_node_t* gnwa_graph_find_max_node(gnwa_node_t** nodes, int32_t n_nodes) {
    // Iterate over all nodes to find the top node
    for (int i = 0; i < n_nodes; ++i) {
        if (nodes[i]->count_prev == 0) {
            return nodes[i];
        }
    }

    fprintf(stderr, "Error:[gnwa] graph has no clear top node!\n");
    exit(1);
}


// Function to destroy a graph
void gnwa_graph_destroy(gnwa_graph_t* graph) {
    for (int i = 0; i < graph->size; ++i) {
        gnwa_node_destroy(graph->nodes[i]);
    }
    graph->max_node = NULL;
    free(graph->nodes);
    graph->nodes = NULL;
    free(graph);
}


// Function to add a node to a graph
uint32_t gnwa_graph_add_node(gnwa_graph_t* graph, gnwa_node_t* node) {
    if (graph->size % 1024 == 0) {
        size_t old_size = graph->size * sizeof(void*);
        size_t increment = 1024 * sizeof(void*);
        if (!(graph->nodes = realloc((void*)graph->nodes, old_size + increment))) {
            fprintf(stderr, "Error:[gnwa] could not allocate memory for graph\n");
            exit(1);
        }
    }
    graph->size++;
    graph->nodes[graph->size - 1] = node;
    return graph->size;
}


// Function to recursively traverse the graph and count the number of paths
uint32_t __gnwa_graph_get_path_count(gnwa_node_t* node) {
    // Check if we have reached the end of the path
    if (node->count_next == 0) {
        return 1;
    }

    uint32_t count = 0;

    // Otherwise follow the path to all of the subsequent nodes
    for (int i = 0; i < node->count_next; ++i) {
        count += __gnwa_graph_get_path_count(node->next[i]);
    }
    return count;
}


// Function to get the count of possible paths
uint32_t gnwa_graph_get_path_count(gnwa_graph_t* graph) {
    // Start at the max node of the graph and take every path recursively
    return __gnwa_graph_get_path_count(graph->max_node);
}


void gnwa_path_build(gnwa_graph_paths_t* paths, gnwa_node_t** nodes, int32_t len) {
    // Check if the end of a path is reached
    if (nodes[len-1]->count_next == 0) {
        gnwa_path_t* path = gnwa_path_create(nodes, len);
        paths->paths[paths->counter] = path;
        paths->counter++;
        return;
    }

    // Iterate over all possible continuing paths
    for (int i = 0; i < nodes[len-1]->count_next; ++i) {
        nodes[len] = nodes[len-1]->next[i];
        gnwa_path_build(paths, nodes, len + 1);
    }
}


// Function to get all paths for a graph
gnwa_graph_paths_t* gnwa_graph_get_paths(gnwa_graph_t* graph) {
    // Allocate memory
    gnwa_graph_paths_t* paths = (gnwa_graph_paths_t*)malloc(sizeof(gnwa_graph_paths_t));
    paths->n_paths = gnwa_graph_get_path_count(graph);
    paths->paths = (gnwa_path_t**)malloc(paths->n_paths * sizeof(gnwa_path_t*));
    paths->counter = 0;

    // Fill the paths into the array
    gnwa_node_t** nodes = (gnwa_node_t**)malloc(graph->size * sizeof(gnwa_node_t*));
    nodes[0] = graph->max_node;
    gnwa_path_build(paths, nodes, 1);
    free(nodes);

    return paths;
}


// Function to destroy the gnwa_graph_paths_t struct (this does not delete the individual paths)
void gnwa_grpah_paths_destroy(gnwa_graph_paths_t* paths) {
    free(paths->paths);
    free(paths);
}


// Function to print a graph in GFA format
void gnwa_graph_print(FILE* file, gnwa_graph_t* graph) {
    // Print the header
    fprintf(file, "H\tVN:Z:1.0\n");

    // Print each node (segment lines)
    for (uint32_t i = 0; i < graph->size; i++) {
        gnwa_node_t* node = graph->nodes[i];
        fprintf(file, "S\t%lu\t%s\n", node->id, node->seq); // Node ID and sequence
    }

    // Print each edge (edge lines)
    for (uint32_t i = 0; i < graph->size; i++) {
        gnwa_node_t* node = graph->nodes[i];

        // Print edges for each next node
        for (int32_t j = 0; j < node->count_next; j++) {
            gnwa_node_t* target_node = node->next[j];
            fprintf(file, "E\t%lu\t%lu\n", node->id, target_node->id); // From node ID to target node ID
        }
    }
}


// Function to allocate memory for an alignment
gnwa_alignment_t* gnwa_alignment_create(const char* read,
                                        gnwa_graph_t* graph,
                                        int32_t score,
                                        gnwa_path_t* path) {
    // Allocate memory
    gnwa_alignment_t* alignment = (gnwa_alignment_t*)malloc(sizeof(gnwa_alignment_t));
    alignment->read = (char*)malloc(strlen(read) + 1);

    // Fill the fields
    strcpy(alignment->read, read);
    alignment->graph = graph;
    alignment->score = score;
    alignment->path = path;

    return alignment;
}


// Functino to destroy an alignment
void gnwa_alignment_destroy(gnwa_alignment_t* alignment) {
    free(alignment->read);
    gnwa_path_destroy(alignment->path);
    gnwa_cigar_destroy(&(alignment->cigar));
    free(alignment);
}


// Function to perform Needleman-Wunsch algorithm
gnwa_alignment_t* gnwa_align(const int8_t* num_read,
                             const int32_t read_len,
                             const int8_t* num_seq,
                             const int32_t seq_len,
                             gnwa_alignment_t* alignment,
                             int8_t* nt_table,
                             int8_t* score_matrix,
                             int8_t gap_open,
                             int8_t gap_extend) {
                                
    // Initialize matrices for scoring and traceback
    int32_t** score_matrix_dp = (int32_t**)malloc((read_len + 1) * sizeof(int32_t*));
    int32_t** trace_matrix = (int32_t**)malloc((read_len + 1) * sizeof(int32_t*));
    for (int i = 0; i <= read_len; ++i) {
        score_matrix_dp[i] = (int32_t*)malloc((seq_len + 1) * sizeof(int32_t));
        // trace_matrix[i] = (int32_t*)malloc((seq_len + 1) * sizeof(int32_t));
        trace_matrix[i] = calloc(1, (seq_len + 1) * sizeof(int32_t));
    }

    // Initialize scoring matrix with gap penalties
    score_matrix_dp[0][0] = 0;
    for (int i = 1; i <= read_len; ++i) {
        score_matrix_dp[i][0] = gap_open + (i - 1) * gap_extend;
        trace_matrix[i][0] = 1; // 1 indicates gap in sequence
    }
    for (int j = 1; j <= seq_len; ++j) {
        score_matrix_dp[0][j] = gap_open + (j - 1) * gap_extend;
        trace_matrix[0][j] = 2; // 2 indicates gap in read
    }

    // Fill the scoring matrix
    for (int i = 1; i <= read_len; ++i) {
        for (int j = 1; j <= seq_len; ++j) {
            int8_t read_base = num_read[i - 1];
            int8_t seq_base = num_seq[j - 1];
            int match_score = score_matrix[read_base * 4 + seq_base]; // Using score matrix

            // Calculate scores for match/mismatch, insertion, and deletion
            int match = score_matrix_dp[i - 1][j - 1] + match_score;
            int delete = score_matrix_dp[i - 1][j] + (trace_matrix[i - 1][j] == 1 ? gap_extend : gap_open);
            int insert = score_matrix_dp[i][j - 1] + (trace_matrix[i][j - 1] == 2 ? gap_extend : gap_open);

            // Choose the maximum score and update the traceback matrix
            score_matrix_dp[i][j] = match;
            trace_matrix[i][j] = 0; // 0 means match/mismatch by default
            if (delete > score_matrix_dp[i][j]) {
                score_matrix_dp[i][j] = delete;
                trace_matrix[i][j] = 1; // 1 means gap in sequence (deletion)
            }
            if (insert > score_matrix_dp[i][j]) {
                score_matrix_dp[i][j] = insert;
                trace_matrix[i][j] = 2; // 2 means gap in read (insertion)
            }
        }
    }

    // Traceback to build the CIGAR string
    int32_t i = read_len;
    int32_t j = seq_len;
    gnwa_cigar_t cigar = {0, NULL}; // Initialize an empty CIGAR struct

    while (i > 0 || j > 0) {
        char op;
        if (trace_matrix[i][j] == 0) {
            // Match/Mismatch
            op = 'M';
            --i;
            --j;
        } else if (trace_matrix[i][j] == 1) {
            // Deletion in sequence
            op = 'D';
            --i;
        } else {
            // Insertion in read
            op = 'I';
            --j;
        }

        // Check if the last CIGAR element has the same type as the current one
        if (cigar.length > 0 && cigar.elements[cigar.length - 1].type == op) {
            // Increment the length of the last element
            cigar.elements[cigar.length - 1].length++;
        } else {
            // Add a new CIGAR element
            cigar.elements = realloc(cigar.elements, (cigar.length + 1) * sizeof(gnwa_cigar_element_t));
            gnwa_cigar_element_t new_element = {op, 1};
            cigar.elements[cigar.length++] = new_element;
        }
    }

    // Reverse the CIGAR elements because traceback builds the alignment from the end to the start
    for (int k = 0; k < cigar.length / 2; ++k) {
        gnwa_cigar_element_t temp = cigar.elements[k];
        cigar.elements[k] = cigar.elements[cigar.length - 1 - k];
        cigar.elements[cigar.length - 1 - k] = temp;
    }

    // Free resources and update alignment score and CIGAR
    alignment->score = score_matrix_dp[read_len][seq_len];
    alignment->cigar = cigar;

    // Free memory used by matrices
    for (int i = 0; i <= read_len; ++i) {
        free(score_matrix_dp[i]);
        free(trace_matrix[i]);
    }
    free(score_matrix_dp);
    free(trace_matrix);

    return alignment;
}


// Function to align a read to a path
gnwa_alignment_t* gnwa_path_align(const char* read,
                                  gnwa_graph_t* graph,
                                  gnwa_path_t* path,
                                  int8_t* nt_table,
                                  int8_t* score_matrix,
                                  uint8_t gap_open,
                                  uint8_t gap_extend) {

    gnwa_alignment_t* alignment = gnwa_alignment_create(read, graph, 0, path);
    int8_t* num_seq = NULL; // Initialize num to NULL
    int32_t seq_len = gnwa_path_get_num_sequence(path, &num_seq); // Pass the address of num
    int32_t read_len = strlen(read); // Get the length of the read
    int8_t* num_read = gnwa_create_num(read, read_len, nt_table); // Create numerical representation of the read


    // Perform the Needleman-Wunsch algorithm on the sequence-read pairs (to be implemented)
    alignment = gnwa_align(num_read,
                            read_len,
                            num_seq,
                            seq_len,
                            alignment,
                            nt_table,
                            score_matrix,
                            -gap_open,
                            -gap_extend);

    // Free allocated memory
    free(num_seq);
    free(num_read);
    return alignment;
}


//  Function to align a read to a graph 
gnwa_alignment_t* gnwa_graph_align(const char* read,
                                    gnwa_graph_t* graph,
                                    int8_t* nt_table,
                                    int8_t* score_matrix,
                                    uint8_t gap_open,
                                    uint8_t gap_extend) {

    // Generate all possible alignment pairs (read + sequence)
    gnwa_graph_paths_t* paths = gnwa_graph_get_paths(graph);

    // Iterate over all paths
    gnwa_alignment_t* alignment = gnwa_path_align(read,
                                                    graph,
                                                    paths->paths[0],
                                                    nt_table,
                                                    score_matrix,
                                                    gap_open,
                                                    gap_extend);
    for (int i = 1; i < paths->n_paths; ++i) {
        gnwa_alignment_t* new_alignment = gnwa_path_align(read,
                                                            graph,
                                                            paths->paths[i],
                                                            nt_table,
                                                            score_matrix,
                                                            gap_open,
                                                            gap_extend);
        // Chosse the alignment with the better score
        if (alignment->score < new_alignment->score) {
            gnwa_alignment_destroy(alignment);
            alignment = new_alignment;
        } else {
            gnwa_alignment_destroy(new_alignment);
        }
    }

    gnwa_grpah_paths_destroy(paths);
    return alignment;
}