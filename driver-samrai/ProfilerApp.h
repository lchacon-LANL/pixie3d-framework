// This class is used to profile functions 
// Copyright © 2010 Mark Berrill. All Rights Reserved. This work is distributed with permission,
// but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
#ifndef included_ProfilerApp
#define included_ProfilerApp


#if defined(WIN32) || defined(_WIN32) || defined(WIN64) || defined(_WIN64)
    #define USE_WINDOWS
#else
    #define USE_LINUX
#endif


#include <iostream>
//using namespace std;

#ifdef USE_WINDOWS
    // Windows
    #include <windows.h>
    #include <string>
    #define TIME_TYPE LARGE_INTEGER
#elif defined(USE_LINUX)
    // Linux
    #include <sys/time.h>
    #include <pthread.h>
    #include <string.h>
    #ifndef LACKS_NAMESPACE
        using namespace std;
    #endif
    #define TIME_TYPE timeval
#else
    #error Unknown OS
#endif


#define BIT_WORD size_t                         // A unsigned integer data type (the larger the word size, the better the performance)
#define TRACE_SIZE 100                          // The maximum number of timers that will be checked for the trace logs
                                                // The actual number of timers is TRACE_SIZE * number of bits of BIT_WORD
#define MAX_TRACE_TRACE 65536                   // The maximum number of stored start and stop times per trace
                                                // Note: this is only used if store_trace is set, and should be a power of 2
                                                // Note: the maximum ammount of memory used per trace is 16*MAX_TRACE_TRACE bytes (plus the trace itself)


/** \class ProfilerApp
  *
  * This class provides some basic timing and profiling capabilities.
  * It works be providing start, stop, and save functions that allow the user to record 
  * the total time required for blocks of code.  The results are written to an ASCII data 
  * file that can be viewed directly or processed.  It is compatible with MPI and SAMRAI.
  * The overhead of the function is ~1us per call for start and stop.  the time for save 
  * depends of the ammount of data to be written and the speed of the writes.  Additional 
  * details can be provided with set_store_trace which records the actual time (from startup)
  * for each start and stop call.  This method significantly increases the memory requirements.
  * Several preprocessor define statement are given for easy incorperation:
  *    PROFILE_START(NAME) - Start profiling a block of code with the given name
  *                          The name must be unique to the file and there must be a corresponding stop statement.
  *    PROFILE_STOP(NAME)  - Stop profiling a block of code with the given name
  *                          The name must match the given start block, and there must only be one PROFILE_STOP for each 
  *                          PROFILE_START.  This records the current line number as the final line number for the block of code.
  *    PROFILE_STOP2(NAME) - Stop profiling a block of code with the given name
  *                          The name must match the given start block.  This is a special stop that does not use the current 
  *                          line as the final line of the block, but is provided for cases where there are multiple exit
  *                          paths for a given function or block of code.
  *    PROFILE_SAVE(FILE)  - Save the results of profiling to a file.  
  * Note that these commands are global and will create a global profiler.  It is possible
  * for a user to create multiple profilers and this should not create any problems, but the 
  * class interface should be used.
  * Example:
  *    #include ProfilerApp.h
  *    void my_function(void *arg) {
  *       PROFILE_START("my function");
  *       int k;
  *       for (int i=0; i<10; i++) {
  *          PROFILE_START("sub call");
  *          // Do some work
  *          if ( special_case1 ) {
  *             PROFILE_STOP2("sub call");
  *             break;
  *          }
  *          if ( special_case2 ) {
  *             PROFILE_STOP2("sub call");
  *             PROFILE_STOP2("my function");
  *             return;
  *          }
  *          // Some more work
  *             PROFILE_STOP2("sub call");
  *       }
  *       PROFILE_STOP("my function");
  *    }
  *    // Some where at the end of the calculation
  *    PROFILE_SAVE("filename");
  */
class ProfilerApp {
public:
    //! Constructor
    ProfilerApp( );

    //! Destructor
    ~ProfilerApp();

   //! Function to start profiling a block of code
    /*!
     * This function starts profiling a block of code until a corresponding stop is called.
     * It is recommended to use PROFILE_START(message) to call this routine.  It will 
     * automatically fill in the file name and the line number.  
     * @param message       Message to uniquely identify the block of code being profiled.
     *                      It must be a unique message to all start called within the same file.
     * @param filename      Name of the file containing the code
     * @param line          Line number containing the start command
     */
    void start( const std::string& message, const std::string& filename, const int line );

    //! Function to stop profiling a block of code
    /*!
     * This function stop profiling a block of code until a corresponding stop is called.
     * It is recommended to use PROFILE_STOP(message) to call this routine.  It will 
     * automatically fill in the file name and the line number.  
     * @param message       Message to uniquely identify the block of code being profiled.
     *                      It must match a start call.
     * @param filename      Name of the file containing the code
     * @param line          Line number containing the stop command
     */
    void stop( const std::string& message, const std::string& filename, const int line );

    //! Function to save the profiling info
    /* Note: .x.timer will automatically be appended to the filename, where x is the rank+1 of the process.
     * Note: .x.trace will automatically be appended to the filename when detailed traces are used.
     * @param filename      File name for saving the results
     */
    void save( const std::string& filename );

    //! Function to change if we are storing detailed trace information (must be called before any start)
    /*  Note: Enabling this option will store the starting and ending time for each call.
     *  This will allow the user to look at the detailed results to get trace information.
     *  However this will significantly increase the memory requirements for any traces
     *  that get called repeatedly and may negitivly impact the performance.
     * @param profile       Do we want to store detailed profiling data
     */
    void set_store_trace(bool profile=false) { if ( id_head==NULL ) { store_trace_data=profile; } }

private:

    // Structure to store the info for a trace log
    struct store_trace {
        int N_calls;                // Number of calls to this block
        BIT_WORD trace[TRACE_SIZE]; // Store the trace
        store_trace *next;          // Pointer to the next entry in the list
        double min_time;            // Store the minimum time spent in the given block (seconds)
        double max_time;            // Store the maximum time spent in the given block (seconds)
        double total_time;          // Store the total time spent in the given block (seconds)
        double *start_time;         // Store when start was called for the given trace (seconds from constructor call)
        double *end_time;           // Store when stop was called for the given trace (seconds from constructor call)
        store_trace() {
            N_calls = 0;
            min_time = 0.0;
            max_time = 0.0;
            total_time = 0.0;
            next = NULL;
            start_time = NULL;
            end_time = NULL;
        }
        // De-constructor used to delete key values
		~store_trace() {
            if ( start_time == NULL )
                delete [] start_time;
            if ( end_time == NULL )
                delete [] end_time;
            start_time = NULL;
            end_time = NULL;
		}
    };

    // Structure to store the timing information for a single block of code
    struct store_timing {
        bool is_active;             // Are we currently running a timer
        std::string message;        // The message to identify the block of code
        std::string filename;       // The file containing the block of code to be timed
        int N_calls;                // Number of calls to this block
        int start_line;             // The starting line for the timer
        int stop_line;              // The ending line for the timer
        BIT_WORD trace[TRACE_SIZE]; // Store the current trace
        double min_time;            // Store the minimum time spent in the given block (seconds)
        double max_time;            // Store the maximum time spent in the given block (seconds)
        double total_time;          // Store the total time spent in the given block (seconds)
        store_trace *trace_head;    // Head of the trace-log list
        store_timing *next;         // Pointer to the next entry in the list
        TIME_TYPE start_time;       // Store when start was called for the given block
        // Constructor used to initialize key values
		store_timing() {
            message = "";
            filename = "";
            start_line = -1;
            stop_line = -1;
			is_active = false;
            min_time = 0.0;
            max_time = 0.0;
            total_time = 0.0;
            N_calls = 0;
            trace_head = NULL;
            next = NULL;
		}
    };
    
    // Structure to store the head of a list for each thread id
    struct store_head {
        #ifdef USE_WINDOWS
            DWORD id;               // The id of the calling thread
        #elif defined(USE_LINUX)
            pthread_t id;           // The id of the calling thread
        #endif
        store_timing *head;         // Pointer to the head of a timing list
        store_head *next;           // Pointer to the next entry in the head list
        store_head *prev;           // Pointer to the previous entry in the head list
        // Constructor used to initialize key values
		store_head() {
            #ifdef USE_WINDOWS
                id = NULL;
            #elif defined(USE_LINUX)
                id = 0;
            #endif
            head = NULL;
            next = NULL;
		}
    };
    
    // Store the profile timers as a linked-list
    store_head *id_head;
    store_head *id_tail;

    // Function to return the appropriate timer block
    store_timing* get_block( const std::string& message, const std::string& filename );

    // Function to return the active timers with the same thread id (it stores the info in binary format)
    void get_active( const int size, BIT_WORD *active );

    // Handle to a mutex lock
    #ifdef USE_WINDOWS
        HANDLE lock;                // Handle to a mutex lock
    #elif defined(USE_LINUX)
        pthread_mutex_t lock;       // Handle to a mutex lock
    #endif
    
    // Misc variables
    bool store_trace_data;          // Do we want to store trace information
    TIME_TYPE construct_time;       // Store when the constructor was called
    TIME_TYPE frequency;            // Clock frequency (only used for windows)
};

extern ProfilerApp global_profiler;
#define PROFILE_START(X)\
    global_profiler.start( X, __FILE__, __LINE__ );
#define PROFILE_STOP(X)\
    global_profiler.stop( X, __FILE__, __LINE__ );
#define PROFILE_START2(X)\
    global_profiler.start( X, __FILE__, -1 );
#define PROFILE_STOP2(X)\
    global_profiler.stop( X, __FILE__, -1 );
#define PROFILE_SAVE(X)\
    global_profiler.save( X );
#define PROFILE_STORE_TRACE(X)\
    global_profiler.set_store_trace( X );

#endif
