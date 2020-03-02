#include <stdio.h> 
#include <stdlib.h> 
#include <time.h>
#include <math.h>

typedef struct AdjListNode { 
    int v; 
    float weight; 
    struct AdjListNode* next; 
} AdjListNode; 
  
typedef struct { 
    AdjListNode* head;
} AdjList ; 

typedef struct { 
    int V; 
    AdjList* array; 
} AdjListGraph; 

AdjListNode* newAdjListNode(int v, float weight) { 
    AdjListNode* newNode = malloc(sizeof(AdjListNode)); 
    newNode->v = v; newNode->weight = weight; newNode->next = NULL; 
    return newNode; 
} 
  
AdjListGraph* createAdjListGraph(int V) { 
    AdjListGraph* G = malloc(sizeof(AdjListGraph)); 
    G->V = V; 
    G->array = malloc(V * sizeof(AdjList)); 
    for (int i = 0; i < V; ++i) 
        G->array[i].head = NULL; 
    return G; 
} 
  
void addEdge(AdjListGraph* G, int src, int dest, float weight) { 
    AdjListNode* newNode = newAdjListNode(dest, weight); 
    newNode->next = G->array[src].head; 
    G->array[src].head = newNode; 
} 
  
typedef struct { 
    int v; 
    float dist; 
} HeapNode; 
  
typedef struct { 
    int d;
    int size;
    int capacity; 
    int* pos;
    HeapNode** array; 
} Heap; 
  
HeapNode* newHeapNode(int v, float dist) { 
    HeapNode* heapNode = malloc(sizeof(HeapNode)); 
    heapNode->v = v; heapNode->dist = dist; 
    return heapNode; 
} 
  
Heap* createHeap(int capacity, int d) { 
    Heap* H = malloc(sizeof(Heap));
    H->array = malloc(capacity * sizeof(HeapNode*)); 
    H->pos = malloc(capacity * sizeof(int)); 
    H->d = d; H->capacity = capacity; H->size = 0; 
    return H; 
} 
  
void swapHeapNode(HeapNode** a, HeapNode** b) { 
    HeapNode* temp = *a; 
    *a = *b; 
    *b = temp; 
} 

int child(int d, int i, int n) {
    return d * i + n;
}
  
void heapify(Heap* H, int i) { 
    int smallest = i; 
    int d = H->d;

    for (int c = child(d, i, 1); c <= child(d, i, d); c++) {
        if (c < H->size && H->array[c]->dist < H->array[smallest]->dist) 
            smallest = c; 
    }
  
    if (smallest != i) { 
  
        H->pos[H->array[smallest]->v] = i; 
        H->pos[H->array[i]->v] = smallest; 
  
        swapHeapNode(&H->array[smallest], &H->array[i]); 
        heapify(H, smallest); 
    } 
} 
  
int isEmpty(Heap* H) { 
    return H->size == 0; 
} 
  
HeapNode* deleteMin(Heap* H) { 
    if (isEmpty(H)) 
        return NULL; 
  
    HeapNode* min = H->array[0]; 
    HeapNode* lastNode = H->array[H->size - 1]; 
    H->array[0] = lastNode; 
  
    H->pos[min->v] = H->size - 1; 
    H->pos[lastNode->v] = 0; 
    --H->size; 

    heapify(H, 0); 
  
    return min; 
} 

int parent(int d, int i) {
    return (i - 1) / d;
}
  
void updateDist(Heap* H, int v, float dist) { 
    int i = H->pos[v]; 
    int d = H->d;
  
    H->array[i]->dist = dist; 
  
    while (i && H->array[i]->dist < H->array[parent(d, i)]->dist) { 
        H->pos[H->array[i]->v] = parent(d, i); 
        H->pos[H->array[parent(d, i)]->v] = i; 

        swapHeapNode(&H->array[i], &H->array[parent(d, i)]); 
        i = parent(d, i); 
    } 
} 
  
int isInHeap(Heap* H, int v) { 
    if (H->pos[v] < H->size) 
        return 1; 
    return 0; 
} 
  
float totalWeight(float arr[], int n) { 
    float weight = 0.;
    for (int i = 1; i < n; ++i) 
        weight += arr[i];
    return weight;
} 

float randf() {
    return (float) rand() / RAND_MAX;
}

float k(int n, int dimension) {
    return 1;
    if (dimension == 4) {
        return pow(n, -0.1);
    }
    return pow(n, -0.5);
}

  
float PrimMST(int V, int dimension) {
    AdjListGraph* G = createAdjListGraph(V);
    float dist[V]; 
    int d = 2;
    if (V > 4) 
        d = k(V, dimension) * (V - 1) / 2;
   
    Heap* H = createHeap(V, d); 

    float threshold = k(V, dimension);

    for (int v = 1; v < V; ++v) { 
        dist[v] = 2.; 
        H->array[v] = newHeapNode(v, dist[v]); 
        H->pos[v] = v; 
    } 
   
    dist[0] = 0; 
    H->array[0] = newHeapNode(0, dist[0]); 
    H->pos[0] = 0; 
    H->size = V; 

    while (!isEmpty(H)) { 
        HeapNode* heapNode = deleteMin(H); 
        int u = heapNode->v; 

        // Add edges only for this vertex.
        for (int i = 0; i < V; i++) {
            if (dimension == 0) {
                float weight = randf();
                if (weight < threshold && i != u) {
                    addEdge(G, u, i, weight); 
                } 
            }
            if (dimension > 0) {
                float sumSquaredDiff = 0;
                for (int i = 0; i < dimension; i++) {
                    float uCoord = randf();
                    float vCoord = randf();
                    sumSquaredDiff += pow(uCoord - vCoord, 2);
                }
                float weight = pow(sumSquaredDiff, .5);
                if (weight < threshold && i != u) {
                    addEdge(G, u, i, weight); 
                } 
            }
        }

        // Explore these edges and update accordingly
        AdjListNode* crawler = G->array[u].head; 
        while (crawler != NULL) { 
            int v = crawler->v; 
            if (isInHeap(H, v) && crawler->weight < dist[v]) { 
                dist[v] = crawler->weight; 
                updateDist(H, v, dist[v]); 
            }
            AdjListNode* temp = crawler;
            crawler = crawler->next;
            free(temp);
        }
        free(heapNode); 
    } 
    free(H->pos); free(H->array); free(H);
    free(G->array); free(G);
    return totalWeight(dist, V);
} 


int main(int argc, char** argv) { 
    srand(time(0));
    int numpoints = atoi(argv[2]);
    int numtrials = atoi(argv[3]);
    int dimension = atoi(argv[4]);
    float sumAverages = 0.;
    for (int i = 0; i < numtrials; i++) {    
        sumAverages += PrimMST(numpoints, dimension);
    }
    float average = sumAverages / numtrials;
    printf("%f %i %i %i\n", average, numpoints, numtrials, dimension);
    return 0;
} 