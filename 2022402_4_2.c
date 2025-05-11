// 2022402 Assignment 4 Task 2
// Parallel Quick Sort using MPI
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <mpi.h>
#include <unistd.h>  // for getopt

// Command line options
struct options {
    int custom_size;  // -n size
    int verbose;      // -v
    int repeat;       // -r count
};

// Function to find median of three values
// 2022402
int medianOfThree(int a, int b, int c) {
    if ((a <= b && b <= c) || (c <= b && b <= a)) return b;
    if ((b <= a && a <= c) || (c <= a && a <= b)) return a;
    return c;
}

// Swap two elements
void swap(int *a, int *b) {
    int t = *a; *a = *b; *b = t;
}

// Partition using median-of-three pivot; returns pivot index
// 2022402
int partition(int arr[], int low, int high) {
    int mid = low + (high - low) / 2;
    int pivot_val = medianOfThree(arr[low], arr[mid], arr[high]);
    
    // Move chosen pivot to end
    if (pivot_val == arr[low])      swap(&arr[low], &arr[high]);
    else if (pivot_val == arr[mid]) swap(&arr[mid], &arr[high]);
    pivot_val = arr[high];

    int i = low - 1;
    for (int j = low; j < high; j++) {
        if (arr[j] <= pivot_val) {
            i++;
            swap(&arr[i], &arr[j]);
        }
    }
    swap(&arr[i+1], &arr[high]);
    return i + 1;
}

// Standard sequential quicksort
void quickSort(int arr[], int low, int high) {
    if (low < high) {
        int pi = partition(arr, low, high);
        quickSort(arr, low, pi - 1);
        quickSort(arr, pi + 1, high);
    }
}

// Merge two sorted arrays
// 2022402
void merge(int arr[], int temp[], int left, int mid, int right) {
    int i = left, j = mid + 1, k = left;
    
    while (i <= mid && j <= right) {
        if (arr[i] <= arr[j])
            temp[k++] = arr[i++];
        else
            temp[k++] = arr[j++];
    }
    
    // Copy remaining elements from first half, if any
    while (i <= mid)
        temp[k++] = arr[i++];
    
    // Copy remaining elements from second half, if any
    while (j <= right)
        temp[k++] = arr[j++];
    
    // Copy back to original array
    for (int x = left; x <= right; x++)
        arr[x] = temp[x];
}

// Verification function to check if array is sorted
int is_sorted(int *a, int n) {
    for (int i = 1; i < n; i++)
        if (a[i] < a[i-1]) return 0;
    return 1;
}

// Helper function to print a small sample of the array
void print_sample(int *arr, int n) {
    printf("First 10 elements: ");
    for (int i = 0; i < (n < 10 ? n : 10); i++) {
        printf("%d ", arr[i]);
    }
    printf("...\n");
}

void print_usage(const char *prog_name) {
    fprintf(stderr, "Usage: %s [options]\n", prog_name);
    fprintf(stderr, "Options:\n");
    fprintf(stderr, "  -n SIZE   Set custom array size\n");
    fprintf(stderr, "  -r COUNT  Number of repeat runs (default: 1)\n");
    fprintf(stderr, "  -v        Verbose mode\n");
    fprintf(stderr, "  -h        Show this help message\n");
}

int main(int argc, char *argv[]) {
    // 2022402
    int rank, size;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // Default sizes to test
    int sizes[] = {1000, 5000, 10000};
    int num_sizes = sizeof(sizes)/sizeof(*sizes);
    
    // Parse command line options (on rank 0)
    struct options opts = {0, 0, 1}; // defaults
    
    if (rank == 0) {
        int opt;
        while ((opt = getopt(argc, argv, "n:r:vh")) != -1) {
            switch (opt) {
                case 'n':
                    opts.custom_size = atoi(optarg);
                    if (opts.custom_size <= 0) {
                        fprintf(stderr, "Array size must be positive\n");
                        MPI_Abort(MPI_COMM_WORLD, 1);
                    }
                    break;
                case 'r':
                    opts.repeat = atoi(optarg);
                    if (opts.repeat <= 0) {
                        fprintf(stderr, "Repeat count must be positive\n");
                        MPI_Abort(MPI_COMM_WORLD, 1);
                    }
                    break;
                case 'v':
                    opts.verbose = 1;
                    break;
                case 'h':
                    print_usage(argv[0]);
                    MPI_Abort(MPI_COMM_WORLD, 0);
                    break;
                default:
                    print_usage(argv[0]);
                    MPI_Abort(MPI_COMM_WORLD, 1);
            }
        }
        
        // Override array size if specified
        if (opts.custom_size > 0) {
            sizes[0] = opts.custom_size;
            num_sizes = 1;
            printf("Using custom array size: %d\n", opts.custom_size);
        }
    }
    
    // Broadcast options to all ranks
    MPI_Bcast(&opts, sizeof(opts), MPI_BYTE, 0, MPI_COMM_WORLD);
    
    // Broadcast the array sizes and count to all ranks
    MPI_Bcast(&num_sizes, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(sizes, num_sizes, MPI_INT, 0, MPI_COMM_WORLD);

    // Run for each problem size
    for (int si = 0; si < num_sizes; si++) {
        // Run multiple times if requested
        for (int run = 0; run < opts.repeat; run++) {
            if (rank == 0 && opts.repeat > 1) {
                printf("\n--- Run %d/%d ---\n", run+1, opts.repeat);
            }
            
            int N = sizes[si];
            
            // Check if array size is sufficient for the number of processes
            if (rank == 0) {
                printf("\n--- N = %d with %d ranks ---\n", N, size);
                if (N < size) {
                    fprintf(stderr,
                        "Error: array size (%d) must be ≥ number of ranks (%d)\n",
                        N, size);
                    MPI_Abort(MPI_COMM_WORLD, 1);
                }
            }

            // Allocate arrays on root only
            int *arr = NULL, *tmp = NULL;
            if (rank == 0) {
                arr = malloc(N * sizeof *arr);
                if (!arr) {
                    perror("malloc for arr");
                    MPI_Abort(MPI_COMM_WORLD, 1);
                }
                
                tmp = malloc(N * sizeof *tmp);
                if (!tmp) {
                    perror("malloc for tmp");
                    free(arr);
                    MPI_Abort(MPI_COMM_WORLD, 1);
                }
                
                // Generate random array
                srand((unsigned)time(NULL) + run); // Different seed for each run
                for (int i = 0; i < N; i++) {
                    arr[i] = rand() % 10000;
                }
                
                if (opts.verbose) {
                    printf("Generated array with %d elements (unsorted)\n", N);
                    print_sample(arr, N);
                }
            }

            // Sequential quicksort timing
            double seq_t = 0;
            if (rank == 0) {
                int *seq_copy = malloc(N * sizeof *seq_copy);
                if (!seq_copy) {
                    perror("malloc for seq_copy");
                    free(arr);
                    free(tmp);
                    MPI_Abort(MPI_COMM_WORLD, 1);
                }
                
                memcpy(seq_copy, arr, N * sizeof *arr);
                
                double t0 = MPI_Wtime();
                quickSort(seq_copy, 0, N - 1);
                seq_t = MPI_Wtime() - t0;
                
                printf("Seq quicksort: %.6f s → %s\n",
                      seq_t, is_sorted(seq_copy, N) ? "OK" : "FAIL");
                
                // Print sample of sorted array for verification
                print_sample(seq_copy, N);
                
                if (opts.verbose) {
                    // Calculate and print statistics about the sequential sorted array
                    long sum = 0;
                    int min = seq_copy[0], max = seq_copy[0];
                    for (int i = 0; i < N; i++) {
                        sum += seq_copy[i];
                        if (seq_copy[i] < min) min = seq_copy[i];
                        if (seq_copy[i] > max) max = seq_copy[i];
                    }
                    double avg = (double)sum / N;
                    printf("Sequential sort stats: min=%d, max=%d, avg=%.2f\n", 
                          min, max, avg);
                }
                
                free(seq_copy);
            }
            
            // Broadcast sequential time to all ranks
            MPI_Bcast(&seq_t, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

            // Calculate local array sizes
            int base = N / size;
            int rem = N % size;
            int local_n = base + (rank < rem ? 1 : 0);
            
            if (opts.verbose && rank == 0) {
                printf("Data distribution:\n");
                for (int r = 0; r < size; r++) {
                    int r_count = base + (r < rem ? 1 : 0);
                    printf("  Rank %d: %d elements (%.1f%%)\n", 
                          r, r_count, (100.0 * r_count) / N);
                }
            }
            
            // Allocate memory for local array
            int *local_arr = malloc(local_n * sizeof *local_arr);
            if (!local_arr) {
                perror("malloc for local_arr");
                MPI_Abort(MPI_COMM_WORLD, 1);
            }
            
            // Initialize arrays for scatter/gather operations
            int *sendcounts = NULL, *displs = NULL;
            if (rank == 0) {
                sendcounts = malloc(size * sizeof *sendcounts);
                displs = malloc(size * sizeof *displs);
                
                if (!sendcounts || !displs) {
                    perror("malloc for sendcounts/displs");
                    MPI_Abort(MPI_COMM_WORLD, 1);
                }
                
                int offset = 0;
                for (int r = 0; r < size; r++) {
                    sendcounts[r] = base + (r < rem ? 1 : 0);
                    displs[r] = offset;
                    offset += sendcounts[r];
                }
            }
            
            // Barrier for clean parallel timing
            MPI_Barrier(MPI_COMM_WORLD);
            double par_start = MPI_Wtime();
            
            // Scatter the data to all processes
            MPI_Scatterv(arr, sendcounts, displs, MPI_INT,
                        local_arr, local_n, MPI_INT,
                        0, MPI_COMM_WORLD);
            
            // Each process sorts its local array
            quickSort(local_arr, 0, local_n - 1);
            
            // Create arrays to hold sorted segments
            int *sorted_arr = NULL;
            if (rank == 0) {
                sorted_arr = malloc(N * sizeof *sorted_arr);
                if (!sorted_arr) {
                    perror("malloc for sorted_arr");
                    MPI_Abort(MPI_COMM_WORLD, 1);
                }
            }
            
            // Gather sorted segments back to root
            MPI_Gatherv(local_arr, local_n, MPI_INT,
                       sorted_arr, sendcounts, displs, MPI_INT,
                       0, MPI_COMM_WORLD);
            
            // Root merges all the segments
            if (rank == 0) {
                // Create segment boundaries for merging
                int *starts = malloc(size * sizeof *starts);
                int *ends = malloc(size * sizeof *ends);
                
                if (!starts || !ends) {
                    perror("malloc for segment boundaries");
                    free(sorted_arr);
                    MPI_Abort(MPI_COMM_WORLD, 1);
                }
                
                // Initialize segment boundaries
                for (int r = 0; r < size; r++) {
                    starts[r] = displs[r];
                    ends[r] = displs[r] + sendcounts[r] - 1;
                }
                
                // Merge segments in a binary tree fashion
                int segments_left = size;
                
                while (segments_left > 1) {
                    int new_segments = 0;
                    
                    for (int i = 0; i < segments_left; i += 2) {
                        if (i + 1 < segments_left) {
                            // Merge two adjacent segments
                            int left = starts[i];
                            int mid = ends[i];
                            int right = ends[i + 1];
                            
                            if (opts.verbose) {
                                printf("Merging segments: [%d-%d] with [%d-%d]\n", 
                                      left, mid, mid+1, right);
                            }
                            
                            merge(sorted_arr, tmp, left, mid, right);
                            
                            // Update for next iteration
                            starts[new_segments] = left;
                            ends[new_segments] = right;
                            new_segments++;
                        } else {
                            // Odd number of segments, keep the last one as is
                            starts[new_segments] = starts[i];
                            ends[new_segments] = ends[i];
                            new_segments++;
                        }
                    }
                    
                    segments_left = new_segments;
                }
                
                // Copy the merged array back to the original array
                memcpy(arr, sorted_arr, N * sizeof *arr);
                
                // Clean up
                free(sorted_arr);
                free(starts);
                free(ends);
            }
            
            // Ensure all processes have completed before stopping the timer
            MPI_Barrier(MPI_COMM_WORLD);
            double par_end = MPI_Wtime();
            double par_t = par_end - par_start;
            
            // Report results
            if (rank == 0) {
                int is_correct = is_sorted(arr, N);
                printf("Par quicksort: %.6f s → %s\n",
                      par_t, is_correct ? "OK" : "FAIL");
                printf("Speedup: %.2fx\n", seq_t / par_t);
                
                // Print sample of sorted array
                print_sample(arr, N);
                
                // Additional stats in verbose mode
                if (opts.verbose) {
                    // Calculate statistics for parallel sorted array
                    long sum = 0;
                    int min = arr[0], max = arr[0];
                    for (int i = 0; i < N; i++) {
                        sum += arr[i];
                        if (arr[i] < min) min = arr[i];
                        if (arr[i] > max) max = arr[i];
                    }
                    double avg = (double)sum / N;
                    printf("Parallel sort stats: min=%d, max=%d, avg=%.2f\n", 
                          min, max, avg);
                    
                    // Performance metrics
                    printf("Performance metrics:\n");
                    printf("  - Elements per rank: %.1f\n", (double)N / size);
                    printf("  - Time per element: %.9f s\n", par_t / N);
                    printf("  - Time * Processors: %.6f processor-seconds\n", par_t * size);
                    printf("  - Parallel efficiency: %.2f%%\n", 100.0 * seq_t / (par_t * size));
                }
                
                // Clean up resources on root
                free(sendcounts);
                free(displs);
                free(arr);
                free(tmp);
            }
            
            // Free local array on all ranks
            free(local_arr);
        }
    }

    MPI_Finalize();
    return 0;
}