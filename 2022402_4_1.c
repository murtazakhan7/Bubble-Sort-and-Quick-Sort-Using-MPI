// 2022402 Assignment 4 Task 1
// Parallel Bubble Sort using MPI
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <mpi.h>
#include <assert.h>
#include <unistd.h>  // for getopt

// Command line options
struct options {
    int custom_size;  // -n size
    int verbose;      // -v
    int repeat;       // -r count
};

// Sequential bubble sort
// 2022402
void bubbleSort(int arr[], int n) {
    for (int i = 0; i < n - 1; i++)
        for (int j = 0; j < n - i - 1; j++)
            if (arr[j] > arr[j + 1]) {
                int tmp = arr[j];
                arr[j] = arr[j + 1];
                arr[j + 1] = tmp;
            }
}

// Merge two sorted arrays
// 2022402
void merge(int arr[], int temp[], int left, int mid, int right) {
    int i = left, j = mid+1, k = left;
    while (i <= mid && j <= right)
        temp[k++] = (arr[i] <= arr[j]) ? arr[i++] : arr[j++];
    while (i <= mid) temp[k++] = arr[i++];
    while (j <= right) temp[k++] = arr[j++];
    for (int x = left; x <= right; x++) arr[x] = temp[x];
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

    // Allow command-line specification of array size
    int array_sizes[] = {1000, 5000, 10000};
    int num_sizes = sizeof(array_sizes)/sizeof(*array_sizes);
    
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
            array_sizes[0] = opts.custom_size;
            num_sizes = 1;
            printf("Using custom array size: %d\n", opts.custom_size);
        }
    }
    
    // Broadcast options to all ranks
    MPI_Bcast(&opts, sizeof(opts), MPI_BYTE, 0, MPI_COMM_WORLD);
    
    // FIX 1: Broadcast the array sizes and count to all ranks
    MPI_Bcast(&num_sizes, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(array_sizes, num_sizes, MPI_INT, 0, MPI_COMM_WORLD); // only broadcast the needed portion

    // Run each problem size the specified number of times
    for (int idx = 0; idx < num_sizes; idx++) {
        for (int run = 0; run < opts.repeat; run++) {
            if (rank == 0 && opts.repeat > 1) {
                printf("\n--- Run %d/%d ---\n", run+1, opts.repeat);
            }
        int N = array_sizes[idx];
        int *arr = NULL, *temp = NULL;
        double seq_t = 0.0, par_start, par_end;

        // 2022402
        if (rank == 0) {
            printf("\n--- N = %d with %d ranks ---\n", N, size);
            if (N < size) {
                fprintf(stderr,
                        "Error: N (%d) < number of ranks (%d)\n",
                        N, size);
                MPI_Abort(MPI_COMM_WORLD, 1);
            }
            arr = malloc(N * sizeof *arr);
            if (!arr) {
                perror("malloc for arr");
                MPI_Abort(MPI_COMM_WORLD, 1);
            }
            
            temp = malloc(N * sizeof *temp);
            if (!temp) {
                perror("malloc for temp");
                free(arr);
                MPI_Abort(MPI_COMM_WORLD, 1);
            }
            
            srand((unsigned)time(NULL));
            for (int i = 0; i < N; i++) arr[i] = rand() % 10000;

            // Time sequential bubble sort
            int *seq = malloc(N * sizeof *seq);
            if (!seq) {
                perror("malloc for seq");
                free(arr);
                free(temp);
                MPI_Abort(MPI_COMM_WORLD, 1);
            }
            
            memcpy(seq, arr, N * sizeof *arr);
            
            // 2022402
            double t0 = MPI_Wtime();
            bubbleSort(seq, N);
            seq_t = MPI_Wtime() - t0;
            
            printf("Seq bubble sort: %.6f s  →  %s\n",
                   seq_t,
                   (is_sorted(seq, N) ? "OK" : "FAIL"));
            
            // Print sample of sorted array for verification
            print_sample(seq, N);
            
            // Additional stats in verbose mode
            if (opts.verbose) {
                // Calculate and print statistics about the array
                long sum = 0;
                int min = seq[0], max = seq[0];
                for (int i = 0; i < N; i++) {
                    sum += seq[i];
                    if (seq[i] < min) min = seq[i];
                    if (seq[i] > max) max = seq[i];
                }
                double avg = (double)sum / N;
                printf("Sequential sort stats: min=%d, max=%d, avg=%.2f\n", 
                       min, max, avg);
            }
            
            free(seq);
        }

        // FIX 2: Broadcast sequential time to all ranks for local speedup calculation
        MPI_Bcast(&seq_t, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

        // Determine local sizes & displacements for data distribution
        // Comment: When size doesn't divide N evenly, the first 'rem' ranks each get one extra element
        int base = N/size, rem = N%size;
        int local_n = base + (rank < rem);
        
        // Report load distribution if verbose
        if (opts.verbose) {
            // Gather all local_n values to rank 0
            int *all_sizes = NULL;
            if (rank == 0) {
                all_sizes = malloc(size * sizeof(int));
                if (!all_sizes) {
                    perror("malloc for all_sizes");
                    MPI_Abort(MPI_COMM_WORLD, 1);
                }
            }
            
            MPI_Gather(&local_n, 1, MPI_INT, all_sizes, 1, MPI_INT, 0, MPI_COMM_WORLD);
            
            if (rank == 0) {
                printf("Load distribution:\n");
                for (int r = 0; r < size; r++) {
                    printf("  Rank %d: %d elements (%.1f%%)\n", 
                           r, all_sizes[r], (100.0 * all_sizes[r]) / N);
                }
                free(all_sizes);
            }
        }
        
        int *local = malloc(local_n * sizeof *local);
        if (!local) {
            perror("malloc for local");
            MPI_Abort(MPI_COMM_WORLD, 1);
        }

        // 2022402
        // FIX 3: Only allocate sendcounts and displs on rank 0
        int *sendcounts = NULL, *displs = NULL;
        if (rank == 0) {
            sendcounts = malloc(size * sizeof *sendcounts);
            if (!sendcounts) {
                perror("malloc for sendcounts");
                MPI_Abort(MPI_COMM_WORLD, 1);
            }
            
            displs = malloc(size * sizeof *displs);
            if (!displs) {
                perror("malloc for displs");
                free(sendcounts);
                MPI_Abort(MPI_COMM_WORLD, 1);
            }
            
            int offset = 0;
            for (int r = 0; r < size; r++) {
                sendcounts[r] = base + (r < rem);
                displs[r] = offset;
                offset += sendcounts[r];
            }
        }
        
        // Verify sendcounts and displs validity on root
        if (rank == 0) {
            // Sanity check to ensure valid count and displacement values
            int total = 0;
            for (int r = 0; r < size; r++) {
                if (sendcounts[r] <= 0) {
                    fprintf(stderr, "Error: sendcounts[%d] = %d is invalid\n", 
                            r, sendcounts[r]);
                    MPI_Abort(MPI_COMM_WORLD, 1);
                }
                if (displs[r] + sendcounts[r] > N) {
                    fprintf(stderr, "Error: displs[%d] + sendcounts[%d] = %d exceeds N = %d\n", 
                            r, r, displs[r] + sendcounts[r], N);
                    MPI_Abort(MPI_COMM_WORLD, 1);
                }
                total += sendcounts[r];
            }
            if (total != N) {
                fprintf(stderr, "Error: Sum of sendcounts = %d doesn't match N = %d\n", 
                        total, N);
                MPI_Abort(MPI_COMM_WORLD, 1);
            }
        }

        // FIX 2: Clean timing of parallel phase using barriers
        // Barrier ensures all processes start parallel section simultaneously
        MPI_Barrier(MPI_COMM_WORLD);
        par_start = MPI_Wtime();

        // Distribute data to all processes
        MPI_Scatterv(arr, sendcounts, displs,
                     MPI_INT, local, local_n,
                     MPI_INT, 0, MPI_COMM_WORLD);
        
        // Each process performs local bubble sort
        bubbleSort(local, local_n);

        // Gather sorted subarrays back to master
        MPI_Gatherv(local, local_n, MPI_INT,
                    arr, sendcounts, displs,
                    MPI_INT, 0, MPI_COMM_WORLD);

        // 2022402
        if (rank == 0) {
            // FIX 4: Corrected merge algorithm for combining sorted segments
            // First create an array to hold the final merged result
            int *merged = malloc(N * sizeof *merged);
            if (!merged) {
                perror("malloc for merged");
                MPI_Abort(MPI_COMM_WORLD, 1);
            }
            
            // Set up arrays for segment boundaries
            int *starts = malloc(size * sizeof *starts);
            int *ends = malloc(size * sizeof *ends);
            if (!starts || !ends) {
                perror("malloc for segment boundaries");
                free(merged);
                if (starts) free(starts);
                MPI_Abort(MPI_COMM_WORLD, 1);
            }
            
            // Initialize segment boundaries
            for (int r = 0; r < size; r++) {
                starts[r] = displs[r];
                ends[r] = displs[r] + sendcounts[r] - 1;
            }
            
            // Merge segments using a modified merge-sort approach
            // This approach merges segments until only one remains
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
                        
                        merge(arr, temp, left, mid, right);
                        
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
            
            // Verify the merged result is fully sorted
            if (!is_sorted(arr, N)) {
                fprintf(stderr, "ERROR: Final merged array is not sorted!\n");
                // Print first 20 elements to help debug
                printf("First 20 elements after merge: ");
                for (int i = 0; i < (N < 20 ? N : 20); i++) {
                    printf("%d ", arr[i]);
                }
                printf("\n");
            }
            
            // Clean up
            free(merged);
            free(starts);
            free(ends);
        }
            
        // FIX 2: Ensure all processes have completed Scatterv/local sort/Gatherv/merge
        // before stopping timer for accurate parallel timing measurement
        MPI_Barrier(MPI_COMM_WORLD);
        par_end = MPI_Wtime();

        if (rank == 0) {
            // Verify & report performance
            double par_t = par_end - par_start;
            printf("Par bubble sort: %.6f s → %s\n",
                   par_t,
                   (is_sorted(arr, N) ? "OK" : "FAIL"));
            printf("Speedup: %.2fx\n", seq_t / par_t);
            
            // Print sample of parallel sorted array for verification
            print_sample(arr, N);
            
            // Additional stats in verbose mode
            if (opts.verbose) {
                // Calculate and print statistics about the parallel sorted array
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
                
                // Timing breakdown
                printf("Performance metrics:\n");
                printf("  - Elements per rank: %.1f\n", (double)N / size);
                printf("  - Time per element: %.9f s\n", par_t / N);
                printf("  - Time * Processors: %.6f processor-seconds\n", par_t * size);
                printf("  - Parallel efficiency: %.2f%%\n", 100.0 * seq_t / (par_t * size));
            }

            // FIX 3: Consistent cleanup for rank 0 
            free(sendcounts);
            free(displs);
            free(arr);
            free(temp);
        }

        // Free local memory on all ranks
        free(local);
        } // End of repeat runs loop
    } // End of array sizes loop

    MPI_Finalize();
    return 0;
}