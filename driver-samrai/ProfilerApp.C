#include "ProfilerApp.h"
#include <stdio.h>

#ifdef USE_WINDOWS
	// Using windows, no MPI
#else
	// Using linux (eventually need to test for SAMRAI)
	//#define USE_SAMRAI
	//#define USE_MPI
#endif

ProfilerApp global_profiler = ProfilerApp();

extern "C" {
    #include "assert.h"
}


#ifdef USE_SAMRAI
    //#include "tbox/Utilities.h"
    //#include "tbox/SAMRAI_MPI.h"
    #include "SAMRAI/tbox/Utilities.h"
    #include "SAMRAI/tbox/SAMRAI_MPI.h"
    #define ERROR_MSG TBOX_ERROR
#elif defined USE_MPI
    #include "mpi.h"
    inline void ERROR_MSG(std::string msg) { return; }
#else
    #define MPI_COMM_WORLD		NULL
    void MPI_Comm_rank(void *comm,int *rank) { *rank = 0; }
    void MPI_Comm_size(void *comm,int *size) { *size = 1; }
    inline void ERROR_MSG(std::string msg) { return; }
#endif


#ifdef USE_WINDOWS
    #define get_time(x) QueryPerformanceCounter(x)
    #define get_diff(start,end,f) (((double)(end.QuadPart-start.QuadPart))/((double)f.QuadPart))
    #define get_frequency(f) QueryPerformanceFrequency(f)
#elif defined(USE_LINUX)
    #define get_time(x) gettimeofday(x,NULL);
    #define get_diff(start,end,f) (((double)end.tv_sec-start.tv_sec)+1e-6*((double)end.tv_usec-start.tv_usec))
    #define get_frequency(f) gettimeofday(f,NULL);
#else
    #error Unknown OS
#endif


#define MAX(a,b) (((a) > (b)) ? (a) : (b))
#define MIN(a,b) (((a) < (b)) ? (a) : (b))


/******************************************************************
* Some inline functions to acquire/release a mutex                *
******************************************************************/
#ifdef USE_WINDOWS
    static inline bool GET_LOCK(HANDLE *lock) {
        int retval = WaitForSingleObject(*lock,INFINITE);
        if ( retval != WAIT_OBJECT_0 ) {
            printf("Error locking mutex\n");
            return true;
        }
	    return false;
    }
    static inline bool RELEASE_LOCK(HANDLE *lock) {
        int retval = ReleaseMutex(*lock);
        if ( retval == 0 ) {
            printf("Error unlocking mutex\n");
            return true;
        }
    	return false;
    }
#else
    static inline bool GET_LOCK(pthread_mutex_t *lock) {
        int retval = pthread_mutex_lock(lock);
        if ( retval == -1 ) {
            printf("Error locking mutex\n");
            return true;
        }
	    return false;
    }
    static inline bool RELEASE_LOCK(pthread_mutex_t *lock) {
        int retval = pthread_mutex_unlock(lock);
        if ( retval == -1 ) {
            printf("Error unlocking mutex\n");
            return true;
        }
	    return false;
    }
#endif


/***********************************************************************
* Consructor                                                           *
***********************************************************************/
ProfilerApp::ProfilerApp() {
    get_frequency( &frequency );
    #ifdef USE_WINDOWS
        lock = CreateMutex (NULL, FALSE, NULL);
    #elif defined(USE_LINUX)
        pthread_mutex_init (&lock,NULL);
    #endif
    id_head = NULL;
    id_tail = NULL;
    store_trace_data = false;
    get_time(&construct_time);
}


/***********************************************************************
* Deconsructor                                                         *
***********************************************************************/
ProfilerApp::~ProfilerApp() {
    store_head *head;
    store_timing *timer, *tmp;
    while ( id_head != NULL ) {
        head = id_head;
        timer = head->head;
        while ( timer != NULL ) {
            store_trace *trace = timer->trace_head;
            while ( trace != NULL ) {
                store_trace *trace_tmp = trace;
                trace = trace->next;
                delete trace_tmp;
            }
            tmp = timer;
            timer = timer->next;
            delete tmp;
        }
        id_head = id_head->next;
        delete head;
    }
}


/***********************************************************************
* Function to start profiling a block of code                          *
***********************************************************************/
void ProfilerApp::start( const std::string& message, const std::string& filename, const int line ) {
    // Get the appropriate timer
    store_timing* timer = get_block( message, filename );
    if ( timer == NULL )
        ERROR_MSG("Failed to get the appropriate timer");
    if ( timer->start_line == -1 ) {
        timer->start_line = line;
    } else if ( line == -1 ) {
        // This is a special case to allow for different start points
    } else if ( timer->start_line != line ) {
        std::string msg = "Multiple start calls with the same message (" + message + ") are not allowed";
        ERROR_MSG(msg);
    }
    if ( timer->is_active ) {
        std::string msg = "Timer is already active, did you forget to call stop? (" + message + ")";
        ERROR_MSG(msg);
    }
    // Get the active trace log
    get_active( TRACE_SIZE, timer->trace );
    // Start the timer 
    timer->is_active = true;
    timer->N_calls++;
    get_time(&timer->start_time);
}


/***********************************************************************
* Function to stop profiling a block of code                           *
***********************************************************************/
void ProfilerApp::stop( const std::string& message, const std::string& filename, const int line ) {
    // Use the current time (minimize the effects of the overhead of the timer)
    TIME_TYPE end_time;
    get_time(&end_time);
    // Get the appropriate timer
    store_timing* timer = get_block( message, filename );
    if ( timer == NULL )
        ERROR_MSG("Failed to get the appropriate timer");
    if ( timer->stop_line == -1 ) {
        timer->stop_line = line;
    } else if ( line == -1 ) {
        // This is a special case to allow for different stop points
    } else if ( timer->stop_line != line ) {
        std::string msg = "Multiple stop calls with the same message (" + message + ") are not allowed";
        ERROR_MSG(msg);
    }
    if ( !timer->is_active ) {
        std::string msg = "Timer is not active, did you forget to call start? (" + message + ")";
        ERROR_MSG(msg);
    }
    timer->is_active = false;
    // Get the active trace log
    BIT_WORD active[TRACE_SIZE];
    get_active( TRACE_SIZE, active );
    // The timer is only a calling timer if it was active before and after the current timer
    if ( memcmp(timer->trace,active,TRACE_SIZE*sizeof(BIT_WORD))!=0 ) {
        for (unsigned int i=0; i<TRACE_SIZE; i++)
            active[i] &= timer->trace[i];
    }
    // Find the trace to save
    store_trace *trace = timer->trace_head;
    while ( trace != NULL) {
        if ( memcmp(trace->trace,active,TRACE_SIZE*sizeof(BIT_WORD))==0 )
            break;
        trace = trace->next;
    }
    if ( trace == NULL ) {
        trace = new store_trace;
        for (unsigned int i=0; i<TRACE_SIZE; i++)
            trace->trace[i] = active[i];
        if ( timer->trace_head == NULL ) {
            timer->trace_head = trace;
        } else {
            store_trace *trace_list = timer->trace_head;
            while ( trace_list->next != NULL)
                trace_list = trace_list->next;
            trace_list->next = trace;
        }
    }
    // Calculate the time elapsed since start was called
    double time = get_diff(timer->start_time,end_time,frequency);
    // Save the starting and ending time if we are storing the detailed traces
    if ( store_trace_data && trace->N_calls<MAX_TRACE_TRACE) {
        // Check if we need to allocate more memory to store the times
        unsigned int size_old, size_new;
        int N = trace->N_calls;
        if ( trace->start_time==NULL ) {
            // We haven't allocated any memory yet
            size_old = 0;
            size_new = 128;
        } else {
            // We want to allocate memory in powers of 2
            size_old = (unsigned int) N;
            // Find the smallest 2^k <= N
            size_new = size_old;
            int N_bits = 0;
            while ( size_new ) {
                N_bits++;
                size_new >>= 1;
            }
            size_new = ((unsigned int)0x01) << N_bits;
            // If size_new==size_old then we have reached the end of the storage and need to allocate more memory
            if ( size_new==size_old )
                size_new *= 2;
            else
                size_new = size_old;
            // Set the minimum ammount of memory to use
            if ( size_new < 128 ) 
                size_new = 128;
            // Stop allocating memory if we reached the limit
            if ( size_new > MAX_TRACE_TRACE ) 
                size_new = MAX_TRACE_TRACE;
        }
        if ( size_old != size_new ) {
            // Expand the trace list
            double *tmp_s = new double[size_new];
            double *tmp_e = new double[size_new];
            for (unsigned int i=0; i<size_old; i++) {
                tmp_s[i] = trace->start_time[i];
                tmp_e[i] = trace->end_time[i];
            }
            if ( trace->start_time!=NULL ) {
                delete [] trace->start_time;
                delete [] trace->end_time;
            }
            trace->start_time = tmp_s;
            trace->end_time = tmp_e;
        }
        // Calculate the time elapsed since the profiler was created
        trace->start_time[N] = time = get_diff(construct_time,timer->start_time,frequency);
        trace->end_time[N] = time = get_diff(construct_time,end_time,frequency);
    }
    // Save the minimum, maximum, and total times
    if ( timer->N_calls == 0 ) {
        timer->min_time = time;
        timer->max_time = time;
    } else {
        timer->max_time = MAX(timer->max_time,time);
        timer->min_time = MIN(timer->min_time,time);
    }
    timer->total_time += time;
    // Save the new time info to the trace
    if ( trace->N_calls == 0 ) {
        trace->min_time = time;
        trace->max_time = time;
    } else {
        trace->max_time = MAX(trace->max_time,time);
        trace->min_time = MIN(trace->min_time,time);
    }
    trace->total_time += time;
    trace->N_calls++;
}



/***********************************************************************
* Function to save the profiling info                                  *
***********************************************************************/
void ProfilerApp::save( const std::string& filename ) {
    // Get the mutex for thread safety (we don't want the list changing while we are saving the data)
    TIME_TYPE end_time;
    get_time(&end_time);
    bool error = GET_LOCK(&lock);
    if ( error )
        return;
    // Open the file(s) for writing
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    char filename_timer[1000], filename_trace[1000];
    sprintf(filename_timer,"%s.%i.timer",filename.c_str(),rank+1);
    sprintf(filename_trace,"%s.%i.trace",filename.c_str(),rank+1);
    FILE *timerFile = fopen(filename_timer,"wb");
    if ( timerFile == NULL ) {
        printf("Error opening file for writing (timer)");
        return;
    }
    FILE *traceFile = NULL;
    if ( store_trace_data ) {
        traceFile = fopen(filename_trace,"wb");
        if ( traceFile == NULL ) {
            printf("Error opening file for writing (trace)");
            return;
        }
    }
    // Create the file header
    fprintf(timerFile,"             Message                 Filename          ID    Start Line  Stop Line  N_calls  Min Time  Max Time  Total Time\n");
    fprintf(timerFile,"---------------------------------------------------------------------------------------------------------------------------\n");
    // Loop through the thread entries, saving the data
    store_head *head = id_head;
    int thread_id = 0;
    while ( head != NULL ) {
        // Loop through the timers saving each one
        store_timing* timer = head->head;
        while ( timer != NULL ) {
            // Get the name of the file (do not include the folder)             
            char filename2[1000];
            int i1 = ((int) timer->filename.find_last_of(47,-1))+1;
            int i2 = ((int) timer->filename.find_last_of(92,-1))+1;
            for (int i=MAX(i1,i2); i<(int) timer->filename.length(); i++)
                filename2[i-MAX(i1,i2)] = timer->filename[i];
            filename2[timer->filename.length()-MAX(i1,i2)] = 0;
            // Get the running times of the timer
            double min_time = timer->min_time;
            double max_time = timer->max_time;
            double tot_time = timer->total_time;
            // If the timer is still running, add the current processing to the totals
            if ( timer->is_active ) {
                double time = get_diff(timer->start_time,end_time,frequency);
                if ( tot_time == 0.0 ) { 
                    min_time = time;
                    max_time = time;
                    tot_time = time;
                } else {
                    min_time = MIN(min_time,time);
                    max_time = MAX(max_time,time);
                    tot_time += time;
                }
            }
            // Save the timer to the file
            fprintf(timerFile,"%24s  %24s",timer->message.c_str(),filename2);
            fprintf(timerFile,"   %4i",thread_id);
            fprintf(timerFile,"   %7i    %7i",timer->start_line,timer->stop_line);
            fprintf(timerFile,"  %8i     %8.3f  %8.3f  %10.3f\n",timer->N_calls,timer->min_time,timer->max_time,timer->total_time);
            timer = timer->next;
        }
        head = head->next;
        thread_id++;
    }
    // Loop through all of the entries, saving the detailed data and the trace logs
    fprintf(timerFile,"\n\n");
    int N_procs;
    MPI_Comm_size(MPI_COMM_WORLD,&N_procs);
    fprintf(timerFile,"<N_procs=%i,id=%i>\n",N_procs,rank);
    head = id_head;
    thread_id = 0;
    get_time(&end_time);
    while ( head != NULL ) {
        store_timing* timer = head->head;
        int id = 0;
        while ( timer != NULL ) {
            // Get the name of the file (do not include the folder)             
            char filename2[1000];
            int i1 = ((int) timer->filename.find_last_of(47,-1))+1;
            int i2 = ((int) timer->filename.find_last_of(92,-1))+1;
            for (int i=MAX(i1,i2); i<(int) timer->filename.length(); i++)
                filename2[i-MAX(i1,i2)] = timer->filename[i];
            filename2[timer->filename.length()-MAX(i1,i2)] = 0;
            // Get the running times of the timer
            double min_time = timer->min_time;
            double max_time = timer->max_time;
            double tot_time = timer->total_time;
            // If the timer is still running, add the current processing time to the totals
            bool add_trace = false;
            BIT_WORD active[TRACE_SIZE];
            double time = 0.0;
            if ( timer->is_active ) {
                time = get_diff(timer->start_time,end_time,frequency);
                min_time = MIN(min_time,time);
                max_time = MAX(min_time,time);
                tot_time += time;
                // Get the active trace log
                add_trace = true;                
                get_active( TRACE_SIZE, active );
                // The timer is only a calling timer if it was active before and after the current timer
                if ( memcmp(timer->trace,active,TRACE_SIZE*sizeof(BIT_WORD))!=0 ) {
                    for (unsigned int i=0; i<TRACE_SIZE; i++)
                        active[i] &= timer->trace[i];
                }
            }
            // Save the timer info
            fprintf(timerFile,"<timer:id=%i,message=%s,file=%s,thread=%i,start=%i,stop=%i,N=%i,min=%e,max=%e,tot=%e>\n",
                id,timer->message.c_str(),filename2,thread_id,timer->start_line,timer->stop_line,
                timer->N_calls,min_time,max_time,tot_time);
            // Store each trace
            store_trace *trace = timer->trace_head;
            while ( trace != NULL ) {
                // Get the running times of the trace
                double trace_min_time = trace->min_time;
                double trace_max_time = trace->max_time;
                double trace_tot_time = trace->total_time;
                // Determine if we need to add the running trace
                if ( add_trace ) {
                    if ( memcmp(trace->trace,active,TRACE_SIZE*sizeof(BIT_WORD))==0 ) {
                        trace_min_time = MIN(trace_min_time,time);
                        trace_max_time = MAX(trace_min_time,time);
                        trace_tot_time += time;
                        add_trace = false;
                    }
                }
                // Save the trace results
                fprintf(timerFile,"<trace:id=%i,N=%i,min=%e,max=%e,tot=%e,active=[ ",
                    id,trace->N_calls,trace_min_time,trace_max_time,trace_tot_time);
                int active_id = 0;
                int BIT_WORD_size = 8*sizeof(BIT_WORD);
                for (unsigned int i=0; i<TRACE_SIZE; i++) {
                    BIT_WORD mask = trace->trace[i];
                    for (int j=0; j<BIT_WORD_size; j++) {
                        if ( (mask&1)!=0 && active_id!= id )
                            fprintf(timerFile,"%i ",active_id);
                        mask >>= 1;
                        active_id++;
                    }
                }
                fprintf(timerFile,"]>\n");
                // Save the detailed trace results (this is a binary file)
                if ( store_trace_data ) { 
                    fprintf(traceFile,"id=%i,N=%i:",id,trace->N_calls);
                    fwrite(trace->start_time,sizeof(double),trace->N_calls,traceFile);
                    fwrite(trace->end_time,sizeof(double),trace->N_calls,traceFile);
                    fprintf(traceFile,"\n");
                }
                // Advance to the next trace
                trace = trace->next;
            }
            // Create a new trace if necessary
            if ( add_trace ) { 
                fprintf(timerFile,"<trace:id=%i,N=%i,min=%e,max=%e,tot=%e,active=[ ",id,1,time,time,time);
                int active_id = 0;
                int BIT_WORD_size = 8*sizeof(BIT_WORD);
                for (unsigned int i=0; i<TRACE_SIZE; i++) {
                    BIT_WORD mask = active[i];
                    for (int j=0; j<BIT_WORD_size; j++) {
                        if ( (mask&1)!=0 && active_id!= id )
                            fprintf(timerFile,"%i ",active_id);
                        mask >>= 1;
                        active_id++;
                    }
                }
                fprintf(timerFile,"]>\n");
                // Save the detailed trace results (this is a binary file)
                if ( store_trace_data ) { 
                    double start_time_trace = time = get_diff(construct_time,timer->start_time,frequency);
                    double end_time_trace = time = get_diff(construct_time,end_time,frequency);
                    fprintf(traceFile,"id=%i,N=%i:",id,1);
                    fwrite(&start_time_trace,sizeof(double),1,traceFile);
                    fwrite(&end_time_trace,sizeof(double),1,traceFile);
                    fprintf(traceFile,"\n");
                }
            }
            timer = timer->next;
            id++;
        }
        head = head->next;
        thread_id++;
    }
    // Close the file(s)
    fclose(timerFile);
    if ( traceFile!=NULL)
        fclose(traceFile);
    // Release the mutex
    RELEASE_LOCK(&lock);
}



/***********************************************************************
* Function to get the timmer for a particular block of code            *
***********************************************************************/
ProfilerApp::store_timing* ProfilerApp::get_block( const std::string& message, const std::string& filename ) {
    // Search for the id (does not need blocking, it either exists or will have this thread's id)
    #ifdef USE_WINDOWS
        DWORD id = GetCurrentThreadId();
    #elif defined(USE_LINUX)
        pthread_t id = pthread_self();
    #endif
    store_head *tmp;
    tmp = id_head;
    while ( tmp != NULL ) {
        if ( id==tmp->id ) {
            // Return the pointer to the appropriate timer
            break;
        }
        tmp = tmp->next;
    }
    // If the id does not exist, create a new id entry (requires blocking)
    store_timing *timer;
    if ( tmp == NULL ) {
        // Create the new id entry
        bool error = GET_LOCK(&lock);
        if ( error )
            return NULL;
        store_head *id_new = new store_head;
        id_new->id = id;
        id_new->head = NULL;
        id_new->next = NULL;
        id_new->prev = NULL;
        if ( id_head == NULL ) {
            id_head = id_new;
            id_tail = id_new;
        } else {
            id_new->prev = id_tail;
            id_tail->next = id_new;
            id_tail = id_new;
        }
        RELEASE_LOCK(&lock);
        // No timer entries, create a new timer
        timer = new store_timing;
        timer->message = message;
        timer->filename = filename;
        timer->start_line = -1;
        timer->stop_line = -1;
        timer->is_active = false;
        timer->N_calls = 0;
        timer->next = NULL;
        id_new->head = timer;
        return timer;
    }
    // Search for the timer (does not need blocking, it either exists or will have this thread's id)
    timer = tmp->head;
    while ( timer != NULL ) {
        if ( message==timer->message && filename==timer->filename ) {
            // Return the pointer to the appropriate timer
            return timer;
        }
        timer = timer->next;
    }
    // Timer not found, create a new timer
    store_timing *prev;
    prev = tmp->head;
    while ( prev->next != NULL )
        prev = prev->next;
    timer = new store_timing;
    timer->message = message;
    timer->filename = filename;
    timer->start_line = -1;
    timer->stop_line = -1;
    timer->is_active = false;
    timer->N_calls = 0;
    timer->next = NULL;
    prev->next = timer;
    return timer;
}


/***********************************************************************
* Function to return the active timers with the same thread id.  It    *
* will store the list of which timers are active in a binary format in *
* an integer array.  Only the timers with the same thread id will be   *
* returned.                                                            *
***********************************************************************/
void ProfilerApp::get_active( const int size, BIT_WORD *active ) {
    // Initialize the active list to 0
    memset(active,0,size*sizeof(BIT_WORD));
    // Search for the id (does not need blocking)
    #ifdef USE_WINDOWS
        DWORD id = GetCurrentThreadId();
    #elif defined(USE_LINUX)
        pthread_t id = pthread_self();
    #endif
    store_head *tmp;
    tmp = id_head;
    while ( tmp != NULL ) {
        if ( id==tmp->id ) {
            // Return the pointer to the appropriate timer
            break;
        }
        tmp = tmp->next;
    }
    if ( tmp == NULL )
        return;
    // Look at all of the timers and see if they are active
    store_timing *timer = tmp->head;
    unsigned int i = 0;
    unsigned int BIT_WORD_size = 8*sizeof(BIT_WORD);
    unsigned int k;
    BIT_WORD mask;
    while ( timer != NULL ) {
        if ( timer->is_active ) {
            // Timer is active, mark it's position
            k = i/BIT_WORD_size;
            mask = 1;
    		mask <<= (i%BIT_WORD_size);
	    	active[k] += mask;
        }
        timer = timer->next;
        i++;
        if ( i >= size*BIT_WORD_size ) {
            // This is a buffer overflow case
            printf("Warning, the active size is not large enough for all the timers, some active timers may be omitted\n");
            break;
        }
    }
}


