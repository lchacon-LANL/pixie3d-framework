#include "ProfilerApp.h"
#include <stdio.h>
#include <iostream>
#include <sstream>


#ifdef USE_WINDOWS
	// Using windows, no MPI
#else
	// Using linux (eventually need to test for SAMRAI)
	//#define USE_SAMRAI
	//#define USE_MPI
#endif

extern ProfilerApp global_profiler = ProfilerApp();

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
    inline void ERROR_MSG(std::string msg) { 
        std::cout << msg << std::endl;
        MPI_Abort(MPI_COMM_WORLD,-1);
        return;
    }
#else
    #define MPI_COMM_WORLD		NULL
    void MPI_Comm_rank(void *comm,int *rank) { *rank = 0; }
    void MPI_Comm_size(void *comm,int *size) { *size = 1; }
    inline void ERROR_MSG(std::string msg) { 
        std::cout << msg << std::endl;
        abort();
        return;
    }
#endif


#ifdef USE_WINDOWS
    #define get_time(x) QueryPerformanceCounter(x)
    #define get_diff(start,end,f) (((double)(end.QuadPart-start.QuadPart))/((double)f.QuadPart))
    #define get_frequency(f) QueryPerformanceFrequency(f)
#elif defined(USE_LINUX)
    #define get_time(x) gettimeofday(x,NULL);
    #define get_diff(start,end,f) (((double)end.tv_sec-start.tv_sec)+1e-6*((double)end.tv_usec-start.tv_usec))
    #define get_frequency(f) (*f=timeval())
#else
    #error Unknown OS
#endif


#define MAX(a,b) (((a) > (b)) ? (a) : (b))
#define MIN(a,b) (((a) < (b)) ? (a) : (b))

template <class type_a, class type_b>
static inline void quicksort2(int n, type_a *arr, type_b *brr);

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
* Inline functions to set or unset the ith bit of the bit array trace  *
***********************************************************************/
static inline void set_trace_bit( unsigned int i, unsigned int N, BIT_WORD *trace ) {
    unsigned int N_bits = 8*sizeof(BIT_WORD);
    unsigned int j = i/N_bits;
    unsigned int k = i%N_bits;
    BIT_WORD mask = ((BIT_WORD)0x1)<<k;
    if ( i < N*N_bits )
        trace[j] |= mask;
}
static inline void unset_trace_bit( unsigned int i, unsigned int N, BIT_WORD *trace ) {
    unsigned int N_bits = 8*sizeof(BIT_WORD);
    unsigned int j = i/N_bits;
    unsigned int k = i%N_bits;
    BIT_WORD mask = ((BIT_WORD)0x1)<<k;
    if ( i < N*N_bits )
        trace[j] &= ~mask;
}


/***********************************************************************
* Inline function to convert the timer id to a string                  *
***********************************************************************/
#define N_BITS_ID 24    // The probability of a collision is ~N^2/2^N_bits (N is the number of timers)
static inline void convert_timer_id( unsigned int id, char* str ) {
    int N_bits = MIN(N_BITS_ID,8*sizeof(unsigned int));
    if ( N_bits <= 9 ) {
        // The id is < 512, store it as a 3-digit number        
        sprintf(str,"%03u",id);
    } else if ( N_bits <= 9 ) {
        // The id is < 2^16, store it as a 4-digit hex
        sprintf(str,"%04x",id);
    } else {
        // We will store the id use the 64 character set { 0-9 a-z A-Z & $ }
        int N = MAX(4,(N_bits+5)/6);    // The number of digits we need to use
        unsigned int tmp1 = id;
        for (int i=N-1; i>=0; i--) {
            unsigned char tmp2 = tmp1%64;
            tmp1 /= 64;
            if ( tmp2 < 10 )
                str[i] = tmp2+48;
            else if ( tmp2 < 36 )
                str[i] = tmp2+(97-10);
            else if ( tmp2 < 62 )
                str[i] = tmp2+(65-36);
            else if ( tmp2 < 63 )
                str[i] = '&';
            else if ( tmp2 < 64 )
                str[i] = '$';
            else
                str[i] = 0;   // We should never use this character
        }
        str[N] = 0;            
    }
}

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
    for (int i=0; i<128; i++)
        thread_head[i] = NULL;
    store_trace_data = false;
    get_time(&construct_time);
    N_threads = 0;
    N_timers = 0;
}


/***********************************************************************
* Deconsructor                                                         *
***********************************************************************/
ProfilerApp::~ProfilerApp() {
    // Delete the thread structures
    for (int i=0; i<128; i++) {
        volatile thread_info *thread = thread_head[i];
        while ( thread != NULL ) {
            // Delete the timers in the thread
            for (int j=0; j<64; j++) {
                store_timer *timer = thread->head[j];
                while ( timer != NULL ) {
                    // Delete the trace logs
                    store_trace *trace = timer->trace_head;
                    while ( trace != NULL ) {
                        store_trace *trace_tmp = trace;
                        trace = trace->next;
                        delete trace_tmp;
                    }
                    store_timer *tmp = timer;
                    timer = timer->next;
                    delete tmp;
                }
            }
            volatile thread_info *thread_next = thread->next;
            delete thread;
            thread = thread_next;
        }
        thread_head[i] = NULL;
    }
    // Delete the global timer info
    for (int i=0; i<128; i++) {
        volatile store_timer_info *timer = timer_table[i];
        while ( timer != NULL ) {
            volatile store_timer_info *timer_next = timer->next;
            delete timer;
            timer = timer_next;
        }
        timer_table[i] = NULL;
    }
}


/***********************************************************************
* Function to start profiling a block of code                          *
***********************************************************************/
void ProfilerApp::start( const std::string& message, const std::string& filename, const int line ) {
    // Get the thread data
    thread_info* thread_data = get_thread_data();
    // Get the appropriate timer
    store_timer* timer = get_block(message,filename,line,-1);
    if ( timer == NULL )
        ERROR_MSG("Failed to get the appropriate timer");
    if ( timer->is_active ) {
        std::stringstream msg;
        msg << "Timer is already active, did you forget to call stop? (" << message << " in " + filename << " at line " << line << ")\n";
        ERROR_MSG(msg.str());
    }
    // Start the timer 
    for (unsigned int i=0; i<TRACE_SIZE; i++)
        timer->trace[i] = thread_data->active[i];
    timer->is_active = true;
    timer->N_calls++;
    set_trace_bit(timer->trace_index,TRACE_SIZE,thread_data->active);
    get_time(&timer->start_time);
}


/***********************************************************************
* Function to stop profiling a block of code                           *
***********************************************************************/
void ProfilerApp::stop( const std::string& message, const std::string& filename, const int line ) {
    // Use the current time (minimize the effects of the overhead of the timer)
    TIME_TYPE end_time;
    get_time(&end_time);
    // Get the thread data
    thread_info* thread_data = get_thread_data();
    // Get the appropriate timer
    store_timer* timer = get_block(message,filename,-1,line);
    if ( timer == NULL )
        ERROR_MSG("Failed to get the appropriate timer");
    if ( !timer->is_active ) {
        std::stringstream msg;
        msg << "Timer is not active, did you forget to call start? (" << message << " in " + filename << " at line " << line << ")\n";
        ERROR_MSG(msg.str());
    }
    timer->is_active = false;
    // Update the active trace log
    unset_trace_bit(timer->trace_index,TRACE_SIZE,thread_data->active );
    // The timer is only a calling timer if it was active before and after the current timer
    BIT_WORD active[TRACE_SIZE];
    for (unsigned int i=0; i<TRACE_SIZE; i++)
        active[i] = thread_data->active[i] & timer->trace[i];
    unsigned int trace_id = get_trace_id( TRACE_SIZE, active );
    // Find the trace to save
    store_trace *trace = timer->trace_head;
    while ( trace != NULL) {
        if ( trace_id==trace->id )
            break;
        trace = trace->next;
    }
    if ( trace == NULL ) {
        trace = new store_trace;
        for (unsigned int i=0; i<TRACE_SIZE; i++)
            trace->trace[i] = active[i];
        trace->id = trace_id;
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
        unsigned int N = trace->N_calls;
        if ( trace->start_time==NULL ) {
            // We haven't allocated any memory yet
            size_old = 0;
            size_new = 1;
        } else {
            // We want to allocate memory in powers of 2
            // The current allocated size is the smallest power of 2 that is >= N
            size_old = 1;
            while ( size_old < N )
                size_old *= 2;
            // Double the storage space (if needed)
            if ( N == size_old )
                size_new = 2*size_old;
            else
                size_new = size_old;
            // Stop allocating memory if we reached the limit
            if ( size_new > MAX_TRACE_TRACE ) 
                size_new = MAX_TRACE_TRACE;
            if ( size_old > MAX_TRACE_TRACE ) 
                size_old = MAX_TRACE_TRACE;
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
        trace->start_time[N] = get_diff(construct_time,timer->start_time,frequency);
        trace->end_time[N] = get_diff(construct_time,end_time,frequency);
    }
    // Save the minimum, maximum, and total times
    if ( timer->N_calls == 1 ) {
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
    // Get the current time in case we need to "stop" and timers
    TIME_TYPE end_time;
    get_time(&end_time);
    // Get the mutex for thread safety (we don't want the list changing while we are saving the data)
    // Note: Because we don't block for most operations in the timer, this is not full proof, but should help
    bool error = GET_LOCK(&lock);
    if ( error )
        return;
    // Get the thread specific data for each thread
    int N_threads2 = N_threads;     // Cache the number of threads since we are holing the lock
    thread_info **thread_data = new thread_info*[N_threads2];
    for (int i=0; i<N_threads2; i++)
        thread_data[i] = NULL;
    for (int i=0; i<128; i++) {
        thread_info *ptr = (thread_info *) thread_head[i];  // It is safe to case to a non-volatile object since we hold the lock
        while ( ptr != NULL ) {
            if ( ptr->thread_num >= N_threads2 )
                ERROR_MSG("Internal error (1)");
            if ( thread_data[ptr->thread_num] != NULL )
                ERROR_MSG("Internal error (2)");
            thread_data[ptr->thread_num] = ptr;
            ptr = (thread_info *) ptr->next;    // It is safe to case to a non-volatile object since we hold the lock
        }
    }
    for (int i=0; i<N_threads2; i++) {
        if ( thread_data[i] == NULL )
            ERROR_MSG("Internal error (3)");
    }
    // Get the timer ids and sort the ids by the total time (maximum value for each thread) to create a global order to save the results
    unsigned int *id_order = new unsigned int[N_timers];
    double *total_time = new double[N_timers];
    for (int i=0; i<N_timers; i++)
        total_time[i] = 0.0;
    int k = 0;
    for (int i=0; i<128; i++) {
        store_timer_info *timer_global = (store_timer_info *) timer_table[i];
        while ( timer_global!=NULL ) {
            id_order[k] = timer_global->id;
            store_timer* timer = NULL;
            for (int thread_id=0; thread_id<N_threads2; thread_id++) {
                thread_info *head = thread_data[thread_id];
                // Search for a timer that matches the current id, and save it
                unsigned int hash = (unsigned int) id_order[k];
                hash *= 0x9E3779B9;     // 2^32*0.5*(sqrt(5)-1)
                unsigned int key = (hash&0xFFFFFFFF)>>26;
                timer = head->head[key];
                while ( timer != NULL ) {
                    if ( timer->id == id_order[k] )
                        break;
                    timer = timer->next;
                }
                if ( timer!=NULL ) {
                    // Get the total running time of the timer
                    total_time[k] = MAX(total_time[k],timer->total_time);
                    // If the timer is still running, add the current processing to the totals
                    if ( timer->is_active ) {
                        double time = get_diff(timer->start_time,end_time,frequency);
                        total_time[k] += time;
                    }
                }
            }
            k++;
            timer_global = (store_timer_info *) timer_global->next;
        }
    }
    if ( k!=N_timers )
        ERROR_MSG("Not all timers were found");
    quicksort2(N_timers,total_time,id_order);
    delete [] total_time;
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
    fprintf(timerFile,"             Message                 Filename        Thread  Start Line  Stop Line  N_calls  Min Time  Max Time  Total Time\n");
    fprintf(timerFile,"---------------------------------------------------------------------------------------------------------------------------\n");
    // Loop through the list of timers, storing the most expensive first
    for (int i=N_timers-1; i>=0; i--) {
        unsigned int id = id_order[i];
        unsigned int hash = ((unsigned int)id)*0x9E3779B9;  // 2^32*0.5*(sqrt(5)-1)
        unsigned int key1 = (hash&0xFFFFFFFF)>>25;   // key1 is 0-127
        unsigned int key2 = (hash&0xFFFFFFFF)>>26;   // key2 is 0-63
        // Search for the global timer info
        store_timer_info *timer_global = (store_timer_info *) timer_table[key1];
        while ( timer_global!=NULL ) {
            if ( timer_global->id == id ) 
                break;
            timer_global = (store_timer_info *) timer_global->next;
        }
        if ( timer_global==NULL )
            ERROR_MSG("Internal error");
        const char* filename2 = timer_global->filename.c_str();
        const char* message = timer_global->message.c_str();
        int start_line = timer_global->start_line;
        int stop_line = timer_global->stop_line;
        // Loop through the thread entries
        for (int thread_id=0; thread_id<N_threads2; thread_id++) {
            thread_info *head = thread_data[thread_id];
            // Search for a timer that matches the current id
            store_timer* timer = head->head[key2];
            while ( timer != NULL ) {
                if ( timer->id == id )
                    break;
                timer = timer->next;
            }
            if ( timer==NULL ) {
                // The current thread does not have a copy of this timer, move on
                continue;
            }
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
            fprintf(timerFile,"%24s  %24s   %4i   %7i    %7i  %8i     %8.3f  %8.3f  %10.3f\n",
                message,filename2,thread_id,start_line,stop_line,timer->N_calls,min_time,max_time,tot_time);
            timer = timer->next;
        }
    }
    // Loop through all of the entries, saving the detailed data and the trace logs
    fprintf(timerFile,"\n\n");
    int N_procs;
    MPI_Comm_size(MPI_COMM_WORLD,&N_procs);
    fprintf(timerFile,"<N_procs=%i,id=%i>\n",N_procs,rank);
    get_time(&end_time);
    char id_str[34];
    // Loop through the list of timers, storing the most expensive first
    for (int i=N_timers-1; i>=0; i--) {
        unsigned int id = id_order[i];
        unsigned int hash = ((unsigned int)id)*0x9E3779B9;  // 2^32*0.5*(sqrt(5)-1)
        unsigned int key1 = (hash&0xFFFFFFFF)>>25;   // key1 is 0-127
        unsigned int key2 = (hash&0xFFFFFFFF)>>26;   // key2 is 0-63
        // Search for the global timer info
        store_timer_info *timer_global = (store_timer_info *) timer_table[key1];
        while ( timer_global!=NULL ) {
            if ( timer_global->id == id ) 
                break;
            timer_global = (store_timer_info *) timer_global->next;
        }
        if ( timer_global==NULL )
            ERROR_MSG("Internal error");
        const char* filename2 = timer_global->filename.c_str();
        const char* message = timer_global->message.c_str();
        int start_line = timer_global->start_line;
        int stop_line = timer_global->stop_line;
        // Loop through the thread entries
        for (int thread_id=0; thread_id<N_threads2; thread_id++) {
            thread_info *head = thread_data[thread_id];
            // Search for a timer that matches the current id
            store_timer* timer = head->head[key2];
            while ( timer != NULL ) {
                if ( timer->id == id )
                    break;
                timer = timer->next;
            }
            if ( timer==NULL ) {
                // The current thread does not have a copy of this timer, move on
                continue;
            }
            // Get the running times of the timer
            double min_time = timer->min_time;
            double max_time = timer->max_time;
            double tot_time = timer->total_time;
            // If the timer is still running, add the current processing time to the totals
            bool add_trace = false;
            double time = 0.0;
            unsigned int trace_id = 0;
            BIT_WORD active[TRACE_SIZE];
            if ( timer->is_active ) {
                add_trace = true;
                time = get_diff(timer->start_time,end_time,frequency);
                min_time = MIN(min_time,time);
                max_time = MAX(min_time,time);
                tot_time += time;
                // The timer is only a calling timer if it was active before and after the current timer
                for (unsigned int i=0; i<TRACE_SIZE; i++)
                    active[i] = head->active[i] & timer->trace[i];
                unset_trace_bit(timer->trace_index,TRACE_SIZE,active);
                trace_id = get_trace_id( TRACE_SIZE, active );
            }
            // Save the timer info
            convert_timer_id(id,id_str);
            fprintf(timerFile,"<timer:id=%s,message=%s,file=%s,thread=%i,start=%i,stop=%i,N=%i,min=%e,max=%e,tot=%e>\n",
                id_str,message,filename2,thread_id,start_line,stop_line,timer->N_calls,min_time,max_time,tot_time);
            // Store each trace
            store_trace *trace = timer->trace_head;
            while ( trace != NULL ) {
                // Get the running times of the trace
                double trace_min_time = trace->min_time;
                double trace_max_time = trace->max_time;
                double trace_tot_time = trace->total_time;
                // Determine if we need to add the running trace
                if ( add_trace ) {
                    if ( trace_id == trace->id ) {
                        trace_min_time = MIN(trace_min_time,time);
                        trace_max_time = MAX(trace_min_time,time);
                        trace_tot_time += time;
                        add_trace = false;
                    }
                }
                // Save the trace results
                convert_timer_id(id,id_str);
                fprintf(timerFile,"<trace:id=%s,thread=%i,N=%i,min=%e,max=%e,tot=%e,active=[ ",
                    id_str,thread_id,trace->N_calls,trace_min_time,trace_max_time,trace_tot_time);
                unsigned int BIT_WORD_size = 8*sizeof(BIT_WORD);
                for (unsigned int i=0; i<TRACE_SIZE; i++) {
                    for (unsigned int j=0; j<BIT_WORD_size; j++) {
                        unsigned int k = i*BIT_WORD_size + j;
                        if ( k == timer->trace_index )
                            continue;
                        BIT_WORD mask = ((BIT_WORD)0x1)<<j;
                        if ( (mask&trace->trace[i])!=0 ) {
                            // The kth timer is active, find the index and write it to the file
                            store_timer* timer_tmp = NULL;
                            for (int m=0; m<64; m++) {
                                timer_tmp = head->head[m];
                                while ( timer_tmp!=NULL ) {
                                    if ( timer_tmp->trace_index==k )
                                        break;
                                    timer_tmp = timer_tmp->next;
                                }
                                if ( timer_tmp!=NULL )
                                    break;
                            }
                            if ( timer_tmp==NULL )
                                ERROR_MSG("Internal Error");
                            convert_timer_id(timer_tmp->id,id_str);
                            fprintf(timerFile,"%s ",id_str);
                        }
                    }
                }
                fprintf(timerFile,"]>\n");
                // Save the detailed trace results (this is a binary file)
                if ( store_trace_data ) { 
                    convert_timer_id(id,id_str);
                    fprintf(traceFile,"id=%s,thread=%i,N=%i:",id_str,thread_id,trace->N_calls);
                    fwrite(trace->start_time,sizeof(double),trace->N_calls,traceFile);
                    fwrite(trace->end_time,sizeof(double),trace->N_calls,traceFile);
                    fprintf(traceFile,"\n");
                }
                // Advance to the next trace
                trace = trace->next;
            }
            // Create a new trace if necessary
            if ( add_trace ) { 
                convert_timer_id(id,id_str);
                fprintf(timerFile,"<trace:id=%s,thread=%i,N=%i,min=%e,max=%e,tot=%e,active=[ ",id_str,thread_id,1,time,time,time);
                unsigned int BIT_WORD_size = 8*sizeof(BIT_WORD);
                for (unsigned int i=0; i<TRACE_SIZE; i++) {
                    for (unsigned int j=0; j<BIT_WORD_size; j++) {
                        unsigned int k = i*BIT_WORD_size + j;
                        if ( k == timer->trace_index )
                            continue;
                        BIT_WORD mask = ((BIT_WORD)0x1)<<j;
                        if ( (mask&active[i])!=0 ) {
                            // The kth timer is active, find the index and write it to the file
                            store_timer* timer_tmp = NULL;
                            for (int m=0; m<64; m++) {
                                timer_tmp = head->head[m];
                                while ( timer_tmp!=NULL ) {
                                    if ( timer_tmp->trace_index==k )
                                        break;
                                    timer_tmp = timer_tmp->next;
                                }
                                if ( timer_tmp!=NULL )
                                    break;
                            }
                            if ( timer_tmp==NULL )
                                ERROR_MSG("Internal Error");
                            convert_timer_id(timer_tmp->id,id_str);
                            fprintf(timerFile,"%s ",id_str);
                        }
                    }
                }
                fprintf(timerFile,"]>\n");
                // Save the detailed trace results (this is a binary file)
                if ( store_trace_data ) { 
                    double start_time_trace = time = get_diff(construct_time,timer->start_time,frequency);
                    double end_time_trace = time = get_diff(construct_time,end_time,frequency);
                    convert_timer_id(id,id_str);
                    fprintf(traceFile,"id=%s,thread=%i,N=%i:",id_str,thread_id,1);
                    fwrite(&start_time_trace,sizeof(double),1,traceFile);
                    fwrite(&end_time_trace,sizeof(double),1,traceFile);
                    fprintf(traceFile,"\n");
                }
            }
        }
    }
    // Close the file(s)
    fclose(timerFile);
    if ( traceFile!=NULL)
        fclose(traceFile);
    // Free temporary memory
    delete [] thread_data;
    delete [] id_order;
    // Release the mutex
    RELEASE_LOCK(&lock);
}


/************************************************************************
* Function to get the data for the current thread                       *
* Note:  If a thread has called this function at some time in the past  *
* then it will be able to return without blocking. When a thread enters *
*  this function for the first time then it will block as necessary.    *
***********************************************************************/
ProfilerApp::thread_info* ProfilerApp::get_thread_data( ) {
    // Get the thread id (as an integer)
    #ifdef USE_WINDOWS
        DWORD thread_id = GetCurrentThreadId();
    #elif defined(USE_LINUX)
        pthread_t thread_id = pthread_self();
    #endif
    // Hash the thread id
    #ifdef USE_WINDOWS
        unsigned int hash = (unsigned int) thread_id;
    #elif defined(USE_LINUX)
	size_t tmp = (size_t) thread_id;
        unsigned int hash = (unsigned int) tmp;
    #endif 
    hash *= 0x9E3779B9;     // 2^32*0.5*(sqrt(5)-1)
    unsigned int key = (hash&0xFFFFFFFF)>>25;
    // Find the first entry with the given key (creating one if necessary)
    if ( thread_head[key]==NULL ) {
        // The entry in the hash table is empty
        // Acquire the lock
        bool error = GET_LOCK(&lock);
        if ( error )
            return NULL;
        // Check if the entry is still NULL
        if ( thread_head[key]==NULL ) {
            // Create a new entry
            thread_head[key] = new thread_info;
            thread_head[key]->id = thread_id;
            thread_head[key]->N_timers = 0;
            thread_head[key]->next = NULL;
            thread_head[key]->thread_num = N_threads;
            N_threads++;
        }
        // Release the lock
        RELEASE_LOCK(&lock);
    }
    volatile thread_info* head = thread_head[key];
    // Find the entry by looking through the list (creating the entry if necessary)
    while ( head->id != thread_id ) {
        // Check if there is another entry to check (and create one if necessary)
        if ( head->next==NULL ) {
            // Acquire the lock
            bool error = GET_LOCK(&lock);
            if ( error )
                return NULL;
            // Check if another thread created an entry while we were waiting for the lock
            if ( head->next==NULL ) {
                // Create a new entry
                thread_info* new_data = new thread_info;
                new_data = new thread_info;
                new_data->id = thread_id;
                new_data->N_timers = 0;
                new_data->next = NULL;
                new_data->thread_num = N_threads;
                N_threads++;
                head->next = new_data;
            }
            // Release the lock
            RELEASE_LOCK(&lock);
        } 
        // Advance to the next entry
        head = head->next;
    }
    // Return the pointer (Note: we no longer need volatile since we are accessing it from the creating thread)
    return (thread_info*) head;
}


/***********************************************************************
* Function to get the timmer for a particular block of code            *
* Note: This function performs some blocking as necessary.             *
***********************************************************************/
ProfilerApp::store_timer* ProfilerApp::get_block( const std::string& message, const std::string& filename1, const int start, const int stop ) {
    // Get the name of the file without the path
    int i1 = ((int) filename1.find_last_of(47))+1;
    int i2 = ((int) filename1.find_last_of(92))+1;
    std::string filename = filename1.substr(MAX(i1,i2),filename1.size());
    // Get the id for the timer
    unsigned int id = get_timer_id(message,filename);
    // Search for the global timer info
    unsigned int hash = (unsigned int) id;
    hash *= 0x9E3779B9;     // 2^32*0.5*(sqrt(5)-1)
    unsigned int key = (hash&0xFFFFFFFF)>>25;   // Get a key 0-127
    if ( timer_table[key]==NULL ) {
        // The global timer does not exist, create it (requires blocking)
        // Acquire the lock
        bool error = GET_LOCK(&lock);
        if ( error )
            return NULL;
        // Check if the entry is still NULL
        if ( timer_table[key]==NULL ) {
            // Create a new entry
            store_timer_info *info_tmp = new store_timer_info;
            info_tmp->id = id;
            info_tmp->start_line = start;
            info_tmp->stop_line = stop;
            info_tmp->message = message;
            info_tmp->filename = filename;
            info_tmp->next = NULL;
            timer_table[key] = info_tmp;
            N_timers++;
        }
        // Release the lock
        RELEASE_LOCK(&lock);
    }
    volatile store_timer_info *info = timer_table[key];
    while ( info->id != id ) {
        // Check if there is another entry to check (and create one if necessary)
        if ( info->next==NULL ) {
            // Acquire the lock
            bool error = GET_LOCK(&lock);
            if ( error )
                return NULL;
            // Check if another thread created an entry while we were waiting for the lock
            if ( info->next==NULL ) {
                // Create a new entry
                store_timer_info *info_tmp = new store_timer_info;
                info_tmp->id = id;
                info_tmp->start_line = start;
                info_tmp->stop_line = stop;
                info_tmp->message = message;
                info_tmp->filename = filename;
                info_tmp->next = NULL;
                info->next = info_tmp;
                N_timers++;
            }
            // Release the lock
            RELEASE_LOCK(&lock);
        } 
        // Advance to the next entry
        info = info->next;
    }
    // Check that the correct timer was found and it is unique
    store_timer_info *info_tmp = (store_timer_info *) info;
    if ( message!=info_tmp->message || filename!=info_tmp->filename ) {
        std::stringstream msg;
        msg << "Error: multiple timers with the same id were detected (" << id << ")\n" << 
            "    " << info_tmp->filename << "   " << info_tmp->message << std::endl << 
            "    " << filename << "   " << message << std::endl;
        ERROR_MSG(msg.str());
    } 
    if ( start==-1 ) {
        // We either are dealing with a stop statement, or the special case for multiple start lines
    } else if ( info->start_line==-1 ) {
        // The timer without a start line, assign it now 
        // Note:  Technically this should be a blocking call, however it is possible to update the start line directly.  
        info->start_line = start;
    } else if ( info->start_line != start ) {
        // Multiple start lines were detected indicating duplicate timers
        std::stringstream msg;
        msg << "Multiple start calls with the same message are not allowed ("
            << message << " in " << filename << " at lines " << start << ", " << info->start_line << ")\n";
        ERROR_MSG(msg.str());
    }
    if ( stop==-1 ) {
        // We either are dealing with a stop statement, or the special case for multiple start lines
    } else if ( info->stop_line==-1 ) {
        // The timer without a start line, assign it now (this requires blocking)
        // Note:  Technically this should be a blocking call, however it is possible to update the stop line directly.  
        info->stop_line = stop;
    } else if ( info->stop_line != stop ) {
        // Multiple start lines were detected indicating duplicate timers
        std::stringstream msg;
        msg << "Multiple start calls with the same message are not allowed ("
            << message << " in " << filename << " at lines " << stop << ", " << info->stop_line << ")\n";
        ERROR_MSG(msg.str());
    }
    // Get the thread-specific data block
    thread_info *thread_data = get_thread_data();
    // Search for the thread-specific timer and create it if necessary (does not need blocking)
    key = (hash&0xFFFFFFFF)>>26;   // Get a key 0-64
    if ( thread_data->head[key]==NULL ) {
        // The timer does not exist, create it
        store_timer *new_timer = new store_timer;
        new_timer->id = id;
        new_timer->is_active = false;
        new_timer->N_calls = 0;
        new_timer->next = NULL;
        new_timer->trace_index = thread_data->N_timers;
        thread_data->N_timers++;
        thread_data->head[key] = new_timer;
    }
    store_timer *timer = thread_data->head[key];
    while ( timer->id != id ) {
        // Check if there is another entry to check (and create one if necessary)
        if ( timer->next==NULL ) {
            store_timer *new_timer = new store_timer;
            new_timer->id = id;
            new_timer->is_active = false;
            new_timer->N_calls = 0;
            new_timer->next = NULL;
            new_timer->trace_index = thread_data->N_timers;
            thread_data->N_timers++;
            timer->next = new_timer;
        } 
        // Advance to the next entry
        timer = timer->next;
    }
    return timer;
}


/***********************************************************************
* Function to return a unique id based on the message and filename.    *
* Note:  We want to return a unique (but deterministic) id for each    *
* filename/message pair.  We want each process or thread to return the *
* same id independent of the other calls.                              *
***********************************************************************/
inline unsigned int ProfilerApp::get_timer_id( const std::string& message, const std::string& filename )
{
    int c;
    // Hash the filename using DJB2
    const char *s = filename.c_str();
    unsigned int hash1 = 5381;
    while((c = *s++)) {
        // hash = hash * 33 ^ c
        hash1 = ((hash1 << 5) + hash1) ^ c;
    }
    // Hash the message using DJB2
    s = message.c_str();
    unsigned int hash2 = 5381;
    while((c = *s++)) {
        // hash = hash * 33 ^ c
        hash2 = ((hash2 << 5) + hash2) ^ c;
    }
    // Combine the two hashes
    unsigned int key = hash1^hash2;
    if ( N_BITS_ID < 8*sizeof(unsigned int) )
        key = (key*0x9E3779B9) >> (8*sizeof(unsigned int)-N_BITS_ID);
    return key;
}


/***********************************************************************
* Function to return a unique id based on the active timer bit array.  *
* This function works by performing a DJB2 hash on the bit array       *
***********************************************************************/
inline unsigned int ProfilerApp::get_trace_id( int N, BIT_WORD *trace ) {
    unsigned int hash = 5381;
    unsigned const char *s = (unsigned char*) trace;
    int N_words = (int)N*sizeof(BIT_WORD)/sizeof(char);
    for (int i=0; i<N_words; i++) {
        int c = (int) s[i];
        // hash = hash * 33 ^ c
        hash = ((hash << 5) + hash) ^ c;
    }
    return hash;
}


/***********************************************************************
* Subroutine to perform a quicksort                                    *
***********************************************************************/
template <class type_a, class type_b>
static inline void quicksort2(int n, type_a *arr, type_b *brr)
{
    bool test;
    int i, ir, j, jstack, k, l, istack[100];
    type_a a, tmp_a;
    type_b b, tmp_b;
    jstack = 0;
    l = 0;
    ir = n-1;
    while (1) {
        if ( ir-l < 7 ) {             // Insertion sort when subarray small enough.
            for ( j=l+1; j<=ir; j++ ) {
                a = arr[j];
                b = brr[j];
                test = true;
                for (i=j-1; i>=0; i--) {
                    if ( arr[i] < a ) {
                        arr[i+1] = a;
                        brr[i+1] = b;
                        test = false;
                        break;
                    }
                    arr[i+1] = arr[i];
                    brr[i+1] = brr[i];
                }
                if ( test ) {
                    i = l-1;
                    arr[i+1] = a;
                    brr[i+1] = b;
                }
            }
            if ( jstack==0 )
                return;
            ir = istack[jstack];    // Pop stack and begin a new round of partitioning.
            l = istack[jstack-1];
            jstack -= 2;
        } else {
            k = (l+ir)/2;           // Choose median of left, center and right elements as partitioning
                                    // element a. Also rearrange so that a(l) ? a(l+1) ? a(ir).
            tmp_a = arr[k];
            arr[k] = arr[l+1];
            arr[l+1] = tmp_a;
            tmp_b = brr[k];
            brr[k] = brr[l+1];
            brr[l+1] = tmp_b;
            if ( arr[l]>arr[ir] ) {
                tmp_a = arr[l];
                arr[l] = arr[ir];
                arr[ir] = tmp_a;
                tmp_b = brr[l];
                brr[l] = brr[ir];
                brr[ir] = tmp_b;
            }
            if ( arr[l+1] > arr[ir] ) {
                tmp_a = arr[l+1];
                arr[l+1] = arr[ir];
                arr[ir] = tmp_a;
                tmp_b = brr[l+1];
                brr[l+1] = brr[ir];
                brr[ir] = tmp_b;
            }
            if ( arr[l] > arr[l+1] ) {
                tmp_a = arr[l];
                arr[l] = arr[l+1];
                arr[l+1] = tmp_a;
                tmp_b = brr[l];
                brr[l] = brr[l+1];
                brr[l+1] = tmp_b;
            }
            // Scan up to find element > a
            j = ir;
            a = arr[l+1];           // Partitioning element.
            b = brr[l+1];
            for (i=l+2; i<=ir; i++) { 
                if ( arr[i]<a ) 
                    continue;
                while ( arr[j]>a )  // Scan down to find element < a.
                    j--;
                if ( j < i )
                    break;          // Pointers crossed. Exit with partitioning complete.
                tmp_a = arr[i];     // Exchange elements of both arrays.
                arr[i] = arr[j];
                arr[j] = tmp_a;
                tmp_b = brr[i];
                brr[i] = brr[j];
                brr[j] = tmp_b;
            }
            arr[l+1] = arr[j];      // Insert partitioning element in both arrays.
            arr[j] = a;
            brr[l+1] = brr[j];
            brr[j] = b;
            jstack += 2;
            // Push pointers to larger subarray on stack, process smaller subarray immediately.
            if ( ir-i+1 >= j-l ) {
                istack[jstack] = ir;
                istack[jstack-1] = i;
                ir = j-1;
            } else {
                istack[jstack] = j-1;
                istack[jstack-1] = l;
                l = i;
            }
        }
    }
}


