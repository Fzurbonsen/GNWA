/*

*/

#include <stdio.h>

#include "gnwa.h"

// Function to read the GFA file and build the graph
gnwa_graph_t* read_gfa_file(const char* filename, int8_t* nt_table, int8_t* score_matrix) {
    FILE* file = fopen(filename, "r");
    if (!file) {
        perror("Could not open GFA file");
        return NULL;
    }

    gnwa_graph_t* graph = gnwa_graph_create(256); // Initialize graph with capacity for 256 nodes
    char line[256]; // Buffer to hold each line of the file

    // Array to store created nodes based on their IDs
    gnwa_node_t* nodes[100] = {NULL}; // Assuming at most 100 nodes, adjust as needed
    uint64_t id;

    while (fgets(line, sizeof(line), file)) {
        // Remove the trailing newline character from the line
        line[strcspn(line, "\n")] = '\0';

        char* token = strtok(line, "\t");
        if (token == NULL) continue;

        // Read segment lines (starting with 'S')
        if (token[0] == 'S') {
            id = strtoull(strtok(NULL, "\t"), NULL, 10); // Read node ID as an integer
            char* sequence = strtok(NULL, "\t");

            // Create a new node
            gnwa_node_t* node = gnwa_node_create(id, sequence, nt_table, score_matrix);
            if (node != NULL) {
                gnwa_graph_add_node(graph, node);
                nodes[id] = node; // Store node in array by its ID
            }
        }
        // Read edge lines (starting with 'E')
        else if (token[0] == 'E') {
            uint64_t from_id = strtoull(strtok(NULL, "\t"), NULL, 10);
            uint64_t to_id = strtoull(strtok(NULL, "\t"), NULL, 10);

            // Add edges between nodes
            if (nodes[from_id] != NULL && nodes[to_id] != NULL) {
                gnwa_node_add_edge(nodes[from_id], nodes[to_id]);
            }
        }
    }

    fclose(file);


    graph->max_node = gnwa_graph_find_max_node(graph->nodes, graph->size);
    return graph;
}



gnwa_graph_t* build_graph(int8_t* nt_table, int8_t* score_matrix) {
    gnwa_graph_t* graph = gnwa_graph_create(4);
    
    gnwa_node_t* node0 = gnwa_node_create(0, "ACTG", nt_table, score_matrix);
    gnwa_node_t* node1 = gnwa_node_create(1, "AGGTC", nt_table, score_matrix);
    gnwa_node_t* node2 = gnwa_node_create(2, "AGGTCCCG", nt_table, score_matrix);
    gnwa_node_t* node3 = gnwa_node_create(3, "AGGTCG", nt_table, score_matrix);

    gnwa_node_add_edge(node0, node1);
    gnwa_node_add_edge(node0, node2);
    gnwa_node_add_edge(node1, node3);
    gnwa_node_add_edge(node2, node3);

    gnwa_graph_add_node(graph, node0);
    gnwa_graph_add_node(graph, node1);
    gnwa_graph_add_node(graph, node2);
    gnwa_graph_add_node(graph, node3);

    return graph;
}

int main() {

    int8_t nt_table[128] = {
		4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
		4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
		4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
		4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
		4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4, 
		4, 4, 4, 4,  3, 3, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
		4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4, 
		4, 4, 4, 4,  3, 3, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4 
	};
    // int8_t score_matrix[16] = {
    //         1,-1,-1,-1,
    //         -1,1,-1,-1,
    //         -1,-1,1,-1,
    //         -1,-1,-1,1
    //         // 1,1,1,1,
    //         // 1,1,1,1,
    //         // 1,1,1,1,
    //         // 1,1,1,1
    // };
    int8_t* score_matrix = gnwa_create_score_matrix(1,1);

    // gnwa_graph_t* graph = build_graph(nt_table, score_matrix);
    gnwa_graph_t* graph = read_gfa_file("tests/test_cases.gfa", nt_table, score_matrix);

    // gnwa_graph_print(stderr, graph);

    gnwa_alignment_t* alignment = gnwa_graph_align("AGGTCCCG",
                                                    graph,
                                                    nt_table,
                                                    score_matrix,
                                                    1, 1);

    fprintf(stderr, "%i\n", alignment->score);
    gnwa_cigar_print(stderr, &(alignment->cigar));
    gnwa_alignment_destroy(alignment);
    free(score_matrix);

    gnwa_graph_destroy(graph);
    fprintf(stderr, "run sucessfull!\n");
    return 0;
}