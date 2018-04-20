// Copyright 2016 UTK JICS AACE
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.


//=============================================================
// 
//  distribute_queries.cpp
//  
//  The function of this code is to read in FASTA formatted file
//  containing the list of sequences to be used as queries in a
//  BLAST search conducted by HPC-BLAST and split the queries
//  the replication groups (MPI level) and thread groups (intra-
//  process, thread level).
//
//
//  Written by - Shane Sawyer
//
//=============================================================


#include <stdlib.h>
#include <string>
#include <cstring>
#include <iostream>
#include <sstream>
#include <fstream>
#include <stdio.h>
#include <limits.h>
#include <stdint.h>
#include <sys/stat.h>
#include <cmath>
#include <time.h>
#include <sys/time.h>
#include <ctype.h>
#include <getopt.h>
//#include <thread>
#include <new>
#include <omp.h>

//#define THRESHOLD 10000
#define MAXBYTES 1073741824
#define max(a,b) ({ __typeof__ (a) _a = (a); \
  __typeof__ (b) _b = (b); \
  _a > _b ? _a : _b; })

#define min(a,b) ({ __typeof__ (a) _a = (a); \
  __typeof__ (b) _b = (b); \
  _a < _b ? _a : _b; })

int compare_length ( const void *a, const void *b );

struct node
{
  uint64_t id;
  uint64_t length;
};

int compare_length ( const void *a, const void *b )
{
  struct node * aa = (struct node*)a;
  struct node * bb = (struct node*)b;

  if ( aa->length > bb->length )
    return -1;
  else if ( aa->length < bb->length )
    return 1;
  else
    return 0;
}

inline bool file_exists(const std::string& name)
{
  struct stat buffer;
  return (stat (name.c_str(), &buffer) == 0);
}

struct membuf : std::streambuf
{
  membuf(char *buf, char *end) { this->setg(buf,buf,end); }
};

void print_usage()
{
  std::cout<<"Usage:\n"<<std::endl;
  std::cout<<"./distribute_queries <options>\n"<<std::endl;
  std::cout<<"     -h     --help               :     Prints the usage message.\n"<<std::endl;
  std::cout<<"     -r     --num_rep_groups     :     REQUIRED - The number of replication groups used in the parallel job."<<std::endl;
  std::cout<<"     -g     --num_thread_groups  :     REQUIRED - The number of thread groups each rank will use."<<std::endl;
  std::cout<<"     -n     --num_threads        :     Optional - The number of threads to use in processing the distribution. Cannot exceed number of replication groups. Defaults to minimum of thread groups and number of cores reported by OpenMP runtime."<<std::endl;
  std::cout<<"     -t     --type               :     REQUIRED - The type of sequences being distributed. This impacts the threshold value; when batches are cut off."<<std::endl;
  std::cout<<"                                                  'a' for amino acids and 'n' for nucleotides."<<std::endl;
  std::cout<<"     -s     --size               :     OPTIONAL - Specify the size of the query batch length threshold. This will override the defaults based on type of sequence."<<std::endl;
  std::cout<<"     -i     --in                 :     REQUIRED - The name of the input query file."<<std::endl;
  std::cout<<"     -o     --out                :     Optional - Name of the output file. Defaults to the input file name with identifications appended for the parallel job."<<std::endl;

  std::cout.flush();
  return;
}

int main(int argc, char* argv[])
{

  // Variable declaration block. ---------------------------------------------------------
  uint64_t num_replication_groups;                    // Number of query distributions per job.
  uint64_t num_thread_groups;                         // Number of query distributions per rank.
  uint64_t i, j;                                      // Loop counters.
  uint64_t num_threads;                               // Number of threads to use in the distribution of queries.
  uint64_t max_num_threads;                           // Max threads supported concurrently in distributing the queries.
  const uint64_t MAX_FILE_NAME_LEN = 276;
  uint64_t total_residues = 0;                        // Total number of residues in the query file.
  uint64_t res_per_rep_grp = 0;                       // Approximate number of residues to assign to each partition.
  bool is_def_line;                                   // Flag to indicate if the '>' character starts a defline (a new query) or appears somewhere in the defline and is not a new query.

  char input_file_name[MAX_FILE_NAME_LEN];
  char output_file_name[MAX_FILE_NAME_LEN];
  uint64_t if_length;                                 // Input file length.
  uint64_t of_length;                                 // Output file length.

  uint64_t num_vbars;                                 // Number of vertical bars encountered. This is for stripping out the query ID tag.
  char qids_unsorted[MAX_FILE_NAME_LEN];              // Name of file containing all the query IDs in the original order.
  std::string print_string;                           // String to print to screen.

  uint64_t length;                                    // Length of the input file in bytes.  
  char *buffer = 0;                                   // Memory buffer for the input file.
  uint64_t num_queries;
  uint64_t *query_offsets = 0;                        // Offsets into the buffer for the start of the queries.
  uint64_t current_query;
  uint64_t *q_lengths = 0;                            // The array that hold the lengths of the queries in the input.
  uint64_t *query_lengths = 0;                        // Pointer to the query lengths.
  uint64_t current_query_length;
  bool start_counting;                                // Flag to indicate length of the sequence may be counted.
  bool look_for_annotation_end;                       // Flag to indicate that we should look for the end of the def line.
  uint64_t partition;                                 // Loop counter over number of partitions applied to the input query file.
  uint64_t *qstart = 0;                               // Array that holds which query the given partition starts with.

  struct node * sequences = 0;                        // Pointer to node structure. Used in sorting the queries by descending length.
  uint64_t ** bins = 0;                               // Bins to put the sorted queries into.
  uint64_t *bin_index = 0;                            // Points to the next open spot.
  uint64_t *current_bin_allocation = 0;               // Stores the number of spaces currently allocated.
  uint64_t realloc_length;                            // How much to add to the current allocation if out of room.
  uint64_t initial_alloc;                             // How much space to initially allocate.
  uint64_t *current_batch_length = 0;                 // Length of the current batch in all bins.
  uint64_t num_search_iterations;                     // Number of iterations the bins in the current partition will take to complete.

  // variables for threaded section
  omp_lock_t token;                                   // Lock variable for restricting access to the standard I/O streams.
  uint64_t tid;                                       // Thread ID.
  uint64_t nthreads;                                  // Number of threads in the thread team spawned by the OpenMP runtime.
  uint64_t td_num_queries;                            // Number of queries the thread owns for a given partition.
  uint64_t qbegin, qend;                              // First and last query in the partition.
  uint64_t td_num_threads;                            // Number of threads as reported in the parallel region.

  uint64_t THRESHOLD;                                 // Number of letters a batch can hold.
  uint64_t user_size;                                 // Number of letters for a batch as specified on the command line.
  // Getopt variables.
  const char* const short_options = "ho:r:g:n:t:s:i:";

  const struct option long_options[] = {
    { "help", 0, NULL, 'h' },
    { "out" , 1, NULL, 'o' },
    { "num_rep_groups", 1, NULL, 'r' },
    { "num_thread_groups", 1, NULL, 'g' },
    { "num_threads", 1, NULL, 'n' },
    { "type", 1, NULL, 't' },
    { "size", 1, NULL, 's' },
    { "in", 1, NULL, 'i' },
    { NULL, 0, NULL, 0 }               };

  int next_option = 0;

  // Flags to indicate reception of arguments.
  bool got_output = false;
  bool got_num_threads = false;
  bool got_input = false;
  bool got_num_groups = false;
  bool got_num_reps = false;
  bool got_type = false;
  bool got_size = false;

  // End variable declaration block. -----------------------------------------------------
  
  // Parse the command line arguments.
  do
    {
      next_option = getopt_long( argc, argv, short_options, long_options, NULL );

      switch( next_option )
	{
	case 'h':
	  print_usage();
	  exit(0);

	case 'o':
	  i = strlen(optarg);
	  if ( i >= (MAX_FILE_NAME_LEN-20) )
	    {
	      std::cerr<<"Error: Output file name length exceeds maximum currently allowed: "<<(MAX_FILE_NAME_LEN-20)<<std::endl;
	      std::cerr.flush();
	      exit(1);
	    }
	  strncpy(output_file_name, optarg, i);
	  output_file_name[i] = '\0';
	  got_output = true;
	  
	  break;

	case 'i':
	  i = strlen(optarg);
	  if ( i >= (MAX_FILE_NAME_LEN-20) )
	    {
	      std::cerr<<"Error: Output file name length exceeds maximum currently allowed: "<<(MAX_FILE_NAME_LEN-20)<<std::endl;
	      std::cerr.flush();
	      exit(1);
	    }
	  strncpy(input_file_name, optarg, i);
	  input_file_name[i]='\0';
	  got_input = true;

	  break;

	case 'r':
	  num_replication_groups = atoi(optarg);

	  got_num_reps = true;

	  break;

	case 'g':
	  num_thread_groups = atoi(optarg);

	  got_num_groups = true;
	  
	  break;

	case 'n':
	  num_threads = atoi(optarg);

	  got_num_threads = true;

	  break;

	case 't':
	  got_type = true;

	  if ( strlen(optarg) != 1 )
	    {
	      std::cerr<<"Error: Unknown sequence type: "<<optarg<<std::endl;
	      std::cerr.flush();
	      exit(1);
	    }

	  if ( optarg[0] != 'a' && optarg[0] != 'n' )
	    {
	      std::cerr<<"Error: Unknown sequence type: "<<optarg<<std::endl;
	      std::cerr.flush();
	      exit(1);
	    }

	  if ( optarg[0] == 'a' )
	    THRESHOLD = 10000;
	  else
	    THRESHOLD = 100000;

	  break;

	case 's':
	  got_size = true;

	  user_size = atoi(optarg);

	  break;
	      
	case '?':
	  std::cerr<<"Error: Invalid command line argument!"<<std::endl;
	  std::cerr.flush();
	  print_usage();
	  exit(1);

	case -1:
	  break;

	default:
	  std::cerr<<"Error: Something unexpected happened while processing the command line options."<<std::endl;
	  std::cerr.flush();
	  exit(1);
	}

    }
  while ( next_option != -1 );

  // Here we set the maximum number of threads for this run.
  max_num_threads = 1;
  
  // OpenMP
  if ( 1 )
    {
      max_num_threads = omp_get_max_threads();
      std::cout<<"OpenMP runtime reports "<<max_num_threads<<" number of threads supported."<<std::endl;
      std::cout.flush();
    }

  // C++ 11 standard.
  //if ( 1 )
  //  {
  //    max_num_threads = std::thread::hardware_concurrency();
  //    std::cout<<"C++11 thread library reports "<<max_num_threads<<" number of threads supported."<<std::endl;
  //    std::cout.flush();
  //  }

  // Perform sanity checks.
  if ( !got_input )
    {
      std::cerr<<"FATAL ERROR: An input file is required."<<std::endl;
      print_usage();
      exit(1);
    }
  else
    {
      if_length = strlen(input_file_name);
    }

  if ( !got_output )
    {
      of_length = if_length;
      strncpy(output_file_name, input_file_name, if_length);
    }
  else
    {
      of_length = strlen(output_file_name);
    }

  if ( got_size )
    {
      // Sanity check. Should be positive.
      if ( user_size < 1 )
	{
	  std::cerr<<"FATAL ERROR: The size of the query batch should be positive. <"<<user_size<<">"<<std::endl;
	  std::cerr.flush();
	  exit(1);
	}

      THRESHOLD = user_size;
      std::cout<<"Using a value of "<<THRESHOLD<<" for the query batch size."<<std::endl;
      std::cout.flush();
    }

  if ( !got_num_reps )
    {
      std::cerr<<"FATAL ERROR: The number of replication groups for the parallel job is required."<<std::endl;
      std::cerr.flush();
      print_usage();
      exit(1);
    }
  else if ( num_replication_groups <= 0 )
    {
      std::cerr<<"Error: The number of replication groups must be positive."<<std::endl;
      std::cerr.flush();
      exit(1);
    }

  if ( !got_num_groups )
    {
      std::cerr<<"FATAL ERROR: The number of thread groups is required. This is the number of query partitions per rank."<<std::endl;
      std::cerr.flush();
      print_usage();
      exit(1);
    }
  else if ( num_thread_groups <= 0 )
    {
      std::cerr<<"Error: The number of thread groups must be positive."<<std::endl;
      std::cerr.flush();
      exit(1);
    }
  else if ( num_thread_groups > 240 )
    {
      std::cerr<<"Error: The number of thread groups exceeds the maximum number of hardware threads for the Xeon Phi 5110p."<<std::endl;
      std::cerr.flush();
      exit(1);
    }

  if ( !got_num_threads )
    {
      // This might be changed so that we always use the most threads possible.

      num_threads = 1;
      std::cout<<"No threads requested for processing the input. Performing the distribution serially."<<std::endl;
    }
  else
    {
      // Got a user value for num_threads to use. Clamp it to the number of replication groups (query partitions).
      max_num_threads = min(max_num_threads,num_replication_groups);

      if ( num_threads <= 0 )
	{
	  std::cerr<<"Error: num_threads must be positive. Setting to 1."<<std::endl;
	  std::cerr.flush();
	  num_threads = 1;
	}
      else if ( num_threads > max_num_threads )
	{
	  std::cerr<<"Error: num_threads greater than the max number of threads <"<<max_num_threads<<">."<<std::endl;
	  std::cerr<<"Resetting thread count to maximum allowed."<<std::endl;
	  std::cerr.flush();
	  num_threads = max_num_threads;
	}
    }

  // Check that the file exists.
  if ( !file_exists(input_file_name) )
    {
      std::cerr<<"FATAL ERROR: "<<input_file_name<<" file was not found!"<<std::endl;
      std::cerr<<"Exiting..."<<std::endl;
      return 1;
    }

  // Open an input file stream.
  std::ifstream ifs(input_file_name , std::ifstream::binary);

  // Get the length of the file.
  ifs.seekg(0,ifs.end);
  length = ifs.tellg();

  std::cout<<"Length of the file is "<<length<<std::endl;

  // Allocate the buffer.
  try
    {
      buffer = new char[length];
    }
  catch(std::bad_alloc &exc)
    {
      std::cerr<<"Allocation for buffer failed. bad_alloc caught."<<std::endl;
      std::cerr.flush();
      exit(1);
    }

  // Copy the file to the buffer.
  ifs.seekg(0,ifs.beg);
  ifs.read(buffer,length);
  ifs.close();

  // Scan the buffer for the number of queries.
  num_queries = 0;
  is_def_line = false;

  // Start the process by looking at the first character of the file.
  if ( buffer[0] == '>' )
    {
      is_def_line = true;
      num_queries++;
    }

  for ( i=1; i < length; ++i )
    {
      // Check to see if we are currently still on a header line (definition line).
      if ( is_def_line && buffer[i] == '\n' )
	is_def_line = false;

      if ( buffer[i] == '>' && not is_def_line )
	{
	  ++num_queries;
	  is_def_line = true;
	}
    }

  // Rescan to find the beginning of each query. Allocate offset array to hold this information.
  try
    {
      query_offsets = new uint64_t[num_queries+1];
    }
  catch(std::bad_alloc &exc)
    {
      std::cerr<<"Allocation for query_offsets failed. bad_alloc caught."<<std::endl;
      std::cerr.flush();
      exit(1);
    }

  current_query = 0;
  is_def_line = false;
  
  // Start the process by looking at the first character of the file.
  if ( buffer[0] == '>' )
    {
      is_def_line = true;
      query_offsets[current_query] = 0;
      ++current_query;
    }

  for ( i=1; i < length; ++i )
    {
      // Check to see if we are currently still on a header line (definition line).
      if ( is_def_line && buffer[i] == '\n' )
	is_def_line = false;

      if ( buffer[i] == '>' && not is_def_line )
	{
	  query_offsets[current_query] = i;
	  ++current_query;
	  is_def_line = true;
	}
    }

  if ( current_query != num_queries )
    {
      std::cerr<<"Error: Found more queries when computing offsets than when counted earlier."<<std::endl;
      std::cerr<<" current_query="<<current_query<<" and num_queries="<<num_queries<<std::endl;
      std::cerr.flush();
      exit(1);
    }
  else
    {
      query_offsets[current_query] = length;
    }

  if ( num_thread_groups > num_queries )
    {
      num_thread_groups = num_queries;
      std::cout<<"Warning: The number of the number of thread groups exceeds the number of queries."<<std::endl;
      std::cout<<"         Resetting the number of thread groups to the number of queries."<<std::endl;
      std::cout.flush();

      if ( num_replication_groups > 1 )
	{
	  std::cout<<"         Use of replication groups is redundant. Setting to 1."<<std::endl;
	  std::cout.flush();
	  num_replication_groups = 1;
	}
    }

  // Compute the query lengths.
  try
    {
      q_lengths = new uint64_t[num_queries+1];
    }
  catch(std::bad_alloc &exc)
    {
      std::cerr<<"Allocation for q_lengths failed. bad_alloc caught."<<std::endl;
      std::cerr.flush();
      exit(1);
    }

  for ( i=0; i <= num_queries; ++i )
    q_lengths[i] = 0;

  start_counting = false;
  look_for_annotation_end = false;
  current_query_length = 0;

  current_query = 0;

  for ( i=0; i < length; ++i )
    {
      if ( buffer[i] == '>' && !look_for_annotation_end )
	{
	  look_for_annotation_end = true;
	  start_counting = false;
	  q_lengths[current_query] = current_query_length;
	  current_query_length = 0;
	  ++current_query;
	}
      else if ( look_for_annotation_end && buffer[i] != '\n' )
	{
	  // Do nothing.
	}
      else if ( buffer[i] == '\n' && look_for_annotation_end )
	{
	  // end of annotation, can count letters in query
	  start_counting = true;
	  look_for_annotation_end = false;
	}
      else if ( start_counting )
	{
	  // Count legitimate characters.
	  if ( isalnum( buffer[i] ) )
	    ++current_query_length;
	}
      else
	{
	  std::cerr<<"Something unexpected happened while counting the lengths of the queries."<<std::endl;
	  std::cerr.flush();
	  break;
	}
    }

  // Set the last length
  q_lengths[current_query] = current_query_length;

  // Set the query_lengths pointer to the first query.
  query_lengths = &( q_lengths[1] );

  // Write out the query identifiers in the original order of the file so that we can use that to reorder the output at the end of the parallel job.
  // Identifiers can be up to 20 characters long.
  
  strncpy(qids_unsorted,"query_ids.unsorted", 18);
  qids_unsorted[18] = '\0';

  std::ofstream ofs;

  ofs.open( qids_unsorted, std::ofstream::out );

  for ( i=0; i < num_queries; ++i )
    {
      // Get the beginning of the ID from the query offset array.
      qbegin = query_offsets[i] + 1; // do not include '>' .

      // Set qend initially to the beginning.
      qend = qbegin;

      while ( (qend-qbegin) < 20 )
	{
	  if ( buffer[qend] == '\n' )
	    break;

	  ++qend;
	}

      if ( (qend-qbegin) < 6 ) // too short
	{
	  std::cerr<<"Error: The query sequence identifier is too short. Query index is "<<i<<std::endl;

	  print_string.assign( &( buffer[qbegin] ), (qend-qbegin) );

	  std::cerr<<"  The def line starts with : "<<print_string<<std::endl;
	  std::cerr.flush();
	}
      if ( (qend-qbegin) > 20 ) // too long
	{
	  std::cerr<<"Error: The query seqeuence identifier is too long. Query index is "<<i<<std::endl;

	  print_string.assign( &( buffer[qbegin] ), (qend-qbegin) );

	  std::cerr<<"  The def line starts with : "<<print_string<<std::endl;
	  std::cerr.flush();
	  qend = qbegin+20;
	}

      // Write out the sequence identifier.
      print_string.assign( &( buffer[qbegin] ), (qend-qbegin) );
      print_string += '\n';

      ofs<<print_string;
    }

  ofs.close();

  // Now we can start the parallel for section on the partitions.

  // Assign the number of queries to each replication group (partition) so that each has roughly the
  // same number of residues.
  qstart = new uint64_t[num_replication_groups+1];
  
  for ( i=0; i < num_queries; ++i )
    {
      total_residues += query_lengths[i];
    }

  res_per_rep_grp = total_residues / num_replication_groups;  // Integer division!
  current_query = 0;

  std::cout<<"Total residues = "<<total_residues<<std::endl;
  std::cout<<"res_per_rep_grp = "<<res_per_rep_grp<<std::endl;

  // Loop over the replication groups up until the last one and assign queries. The last one gets the remainder.
  
  // Degenerate case:
  if ( num_replication_groups == 1 )
    qstart[0]=0; // the only rep grp gets all the queries.
  
  for ( i=0; i < (num_replication_groups-1); ++i )
    {
      uint64_t current_residues = 0;

      qstart[i] = current_query;

      while ( current_residues < res_per_rep_grp )
	{
	  current_residues += query_lengths[current_query];
	  ++current_query;
	  
	  if ( current_query >= num_queries )
	    {
	      // we have exceeded the available queries and must stop before we go past our buffer.
	      std::cerr<<"A problem was encountered in assigning queries to replication groups. On replication group "
		       <<i<<"."<<std::endl;
	      std::cerr.flush();
	      exit(1);
	    }
	}
      std::cout<<"Replication group "<<i<<" has "<<current_residues<<std::endl;
    }

  if ( num_replication_groups > 1 )
    qstart[num_replication_groups-1] = current_query;

  qstart[num_replication_groups] = num_queries;

  for ( i=0; i <= num_replication_groups; i++ )
    {
      std::cout<<"qstart["<<i<<"] = "<<qstart[i]<<std::endl;
      std::cout.flush();
    }
  std::cout.flush();

  omp_set_num_threads(num_threads);

  omp_init_lock(&token);

  #pragma omp parallel for default(none) shared(length,buffer,num_queries,query_offsets,query_lengths,num_replication_groups, \
						num_thread_groups,token,std::cout,std::cerr,of_length,output_file_name,THRESHOLD, \
						qstart ) \
                                         private(i,j,tid,td_num_queries,qbegin,qend,sequences,bins, \
						 bin_index,current_bin_allocation,realloc_length, \
						 initial_alloc,current_batch_length,current_query, \
						 num_search_iterations,current_query_length, \
						 start_counting,look_for_annotation_end,td_num_threads)
  for ( partition=0; partition < num_replication_groups; ++partition )
    {

      // Get thread information.
      td_num_threads = omp_get_num_threads();
      tid = omp_get_thread_num();

      // Determine how many queries this partition has. Then find the beginning and ending query inidices for the partition.
      //td_num_queries = uint64_t( ceil(double(num_queries) / double(num_replication_groups)) );

      //qbegin = partition * td_num_queries;
      //qend = (partition+1)*td_num_queries;

      qbegin = qstart[partition];
      qend   = qstart[partition+1];
      td_num_queries = qend-qbegin;
      /*
      if ( partition == (num_replication_groups-1) )
	{
	  td_num_queries = num_queries - ( td_num_queries*(num_replication_groups-1) );
	  qend = num_queries;
	}
      */
      // Initialize the node structure array for sorting.
      try
	{
	  sequences = new node[td_num_queries];
	}
      catch(std::bad_alloc &exc)
	{
	  omp_set_lock(&token);
	  std::cerr<<"Allocation for the sequences array failed for thread "<<tid<<" working on partition "<<partition;
	  std::cerr<<". bad_alloc caught."<<std::endl;
	  std::cerr.flush();
	  omp_unset_lock(&token);
	  exit(1);
	}

      for ( i=0; i < td_num_queries; ++i )
	{
	  sequences[i].id = (i+qbegin);
	  sequences[i].length = query_lengths[(i+qbegin)];
	}

      // Sort the nodes by sequence length.
      //std::cout<<"Sorting the queries by descending length..."<<std::endl;

      qsort ( sequences , td_num_queries, sizeof(struct node), compare_length );

      //std::cout<<"Sorted queries."<<std::endl;

      // Allocate the bins. All bins allocation will be done with malloc so realloc can be used.
      bins = (uint64_t **) malloc( num_thread_groups * sizeof(uint64_t*) );

      if ( bins == 0 )
	{
	  omp_set_lock(&token);
	  std::cerr<<"Allocation for the bins array failed for thread "<<tid<<" working on partition "<<partition;
	  std::cerr<<". bad_alloc caught."<<std::endl;
	  std::cerr.flush();
	  omp_unset_lock(&token);
	  exit(1);
	}

      // Compute the initial allocation.
      initial_alloc = ( (td_num_queries/num_thread_groups) / 2 );

      // Compute the reallocation increment.
      realloc_length = max( 100 , (initial_alloc/10) );

      for ( i=0; i < num_thread_groups; ++i )
	{     
	  bins[i] = (uint64_t *) malloc( initial_alloc * sizeof(uint64_t) );
	  
	  if ( bins[i] == 0 )
	    {
	      omp_set_lock(&token);
	      std::cerr<<"Allocation for bins["<<i<<"] array failed for thread "<<tid<<" working on partition "<<partition;
              std::cerr<<". bad_alloc caught."<<std::endl;
	      std::cerr.flush();
	      omp_unset_lock(&token);
	      exit(1);
	    }
	}
      
      // Allocate bin trackers -- last added element, current allocation amount.
      bin_index = (uint64_t*) malloc( num_thread_groups * sizeof(uint64_t) );
      
      if ( bin_index == 0 )
	{
	  omp_set_lock(&token);
	  std::cerr<<"Allocation for the bin_index array failed for thread "<<tid<<" working on partition "<<partition;
	  std::cerr<<". bad_alloc caught."<<std::endl;
	  std::cerr.flush();
	  omp_unset_lock(&token);
	  exit(1);
	} 

      for ( i=0; i < num_thread_groups; ++i )
	bin_index[i] = 0;

      current_bin_allocation = (uint64_t*) malloc( num_thread_groups * sizeof(uint64_t) );

      if ( current_bin_allocation == 0 )
	{
	  omp_set_lock(&token);
	  std::cerr<<"Allocation for the current_bin_allocation array failed for thread "<<tid<<" working on partition "<<partition;
	  std::cerr<<". bad_alloc caught."<<std::endl;
	  std::cerr.flush();
	  omp_unset_lock(&token);
	  exit(1);
	}
      
      for ( i=0; i < num_thread_groups; ++i )
	current_bin_allocation[i] = initial_alloc;

      current_batch_length = (uint64_t*) malloc( num_thread_groups * sizeof(uint64_t) );

      if ( current_batch_length == 0 )
	{
	  omp_set_lock(&token);
	  std::cerr<<"Allocation for the current_batch_length array failed for thread "<<tid<<" working on partition "<<partition;
	  std::cerr<<". bad_alloc caught."<<std::endl;
	  std::cerr.flush();
	  omp_unset_lock(&token);
	  exit(1);
	}
      
      for ( i=0; i < num_thread_groups; ++i )
	current_batch_length[i] = 0;

      num_search_iterations = 0;

      current_query = 0;

      bool *bin_state = 0;
      try
	{
	  bin_state = new bool[num_thread_groups];
	}
      catch(std::bad_alloc &exc)
	{
	  omp_set_lock(&token);
	  std::cerr<<"Allocation for bin_state failed for thread "<<tid<<" working on partition "<<partition;
	  std::cerr<<". bad_alloc caught."<<std::endl;
	  std::cerr.flush();
	  omp_unset_lock(&token);
	  exit(1);
	}

      uint64_t *total_batch_length = 0;
      try
	{
	  total_batch_length = new uint64_t[num_thread_groups];
	}
      catch(std::bad_alloc &exc)
	{
	  omp_set_lock(&token);
	  std::cerr<<"Allocation for the sequences array failed for thread "<<tid<<" working on partition "<<partition;
	  std::cerr<<". bad_alloc caught."<<std::endl;
	  std::cerr.flush();
	  omp_unset_lock(&token);
	  exit(1);
	}

      for ( i=0; i < num_thread_groups; ++i )
	{
	  total_batch_length[i] = 0;
	}

      uint64_t td_num_letters = 0;
      
      for ( i=0; i < num_thread_groups; ++i )
	{
	  bin_state[i] = true;
	}
      // true -- room for more
      // false -- exceeded the limit for this batch group

      bool work_to_do = true;
      bool assign_query = true;
      bool reset_state = false;
      uint64_t least_full_bin = 0;
      
      for ( i=0; i < td_num_queries; ++i )
	{
	  assign_query = true;

	  if ( reset_state )
	    {
	      for ( j=0; j < num_thread_groups; ++j )
		bin_state[j] = true;

	      for ( j=0; j < num_thread_groups; ++j )
		current_batch_length[j] = 0;

	      ++num_search_iterations;
	    }

	  // who has the least sequence letters?
      
	  // Assign to the first bin that is still accepting.
	  for ( j=0; j < num_thread_groups; ++j )
	    {
	      if ( bin_state[j] )
		{
		  least_full_bin = j;
		  break;
		}
	    }
      
	  // sanity check -- if we didn't find a valid bin with a minimum, we have to reset the state.
	  if ( j == (num_thread_groups+1) )
	    {
	      --i; // we haven't assigned this query so try again on the next iteration.
	      reset_state = true;
	      continue;
	    }
      
	  // find the actual least bin still accecpting new queries.
	  for ( j=(least_full_bin+1); j < num_thread_groups; ++j )
	    {
	      if ( ( current_batch_length[j] < current_batch_length[least_full_bin] ) && bin_state[j] )
		least_full_bin = j;
	    }

	  // Now assign the query to this bin. Reallocate the bin if too small.

	  if ( (bin_index[least_full_bin]+1) >= current_bin_allocation[least_full_bin] ) // Realloc
	    {
	      bins[least_full_bin] = (uint64_t*) realloc( (void*)(bins[least_full_bin]), ( current_bin_allocation[least_full_bin] +
										      realloc_length )*sizeof(uint64_t) );

	      if ( bins[least_full_bin] == 0 ) // Realloc failed.
		{
		  omp_set_lock(&token);
		  std::cerr<<"Reallocation for bins["<<least_full_bin<<"] array failed for thread "<<tid<<" working on partition ";
		  std::cerr<<partition<<". bad_alloc caught."<<std::endl;
		  std::cerr.flush();
		  omp_unset_lock(&token);
		  exit(1);
		}

	      current_bin_allocation[least_full_bin] += realloc_length;
	    }

	  bins[least_full_bin][bin_index[least_full_bin]] = sequences[i].id;
	  bin_index[least_full_bin] += 1;
	  current_batch_length[least_full_bin] += sequences[i].length ;
	  total_batch_length[least_full_bin] += sequences[i].length;
	  td_num_letters += sequences[i].length;

	  if ( current_batch_length[least_full_bin] > (THRESHOLD - 1) )
	    bin_state[least_full_bin] = false;
	  
	  // if all bins are full, trigger them to be reset
	  reset_state = true;
	  for ( j=0; j < num_thread_groups; ++j )
	    {
	      if ( bin_state[j] ) // if the any of the bin states are true, at least one bin is still accepting more queries.
		reset_state = false;
	    }
	}

      // Print out the statistics for this partition.
      omp_set_lock(&token);
      std::cout<<"-------------------------Statistics for partition "<<partition<<"-------------------------"<<std::endl;
      std::cout<<"Processed by thread "<<tid<<std::endl;
      std::cout<<"  Number of queries in partition - "<<td_num_queries<<std::endl;
      std::cout<<"  Number of letters in partition - "<<td_num_letters<<std::endl;
      for ( i=0; i < num_thread_groups; ++i )
	std::cout<<"    Thread group/bin has "<<(total_batch_length[i])<<" letters."<<std::endl;
      std::cout<<"  The partition will require "<<(num_search_iterations+1)<<" search iterations in the main BLAST loop."<<std::endl;
      std::cout<<std::endl;
      std::cout.flush();
      omp_unset_lock(&token);

      // Now we can write out the partition to disk. This will consist of a sequence of files; one for each bin (thread group).
      char td_out_file_name[MAX_FILE_NAME_LEN];
      std::ofstream ofs;
      uint64_t bytes_written = 0;

      // Now write out each bin for every thread group in the current partition.
      for ( i=0; i < num_thread_groups; ++i )
	{
	  // Create the file name.
	  strncpy(td_out_file_name, output_file_name, of_length);
	  sprintf( &(td_out_file_name[of_length]), ".%d.%d",int(partition),int(i));

	  // Open the stream for writing.
	  ofs.open( td_out_file_name, std::ofstream::out | std::ofstream::binary );

	  uint64_t qlen = 0;            // the length of the query. Based on the query_offsets array.
	  uint64_t qidx = 0;            // the index of the query in the global arrays.

	  for ( j=0; j < bin_index[i]; ++j )
	    {
	      qidx = bins[i][j];

	      qlen = query_offsets[qidx+1] - query_offsets[qidx];

	      ofs.write( &(buffer[ query_offsets[qidx] ]), qlen );

	      bytes_written += qlen;

	      if ( bytes_written > MAXBYTES )
		{
		  omp_set_lock(&token);
		  std::cerr<<"ERROR: In partition "<<partition<<", the file for thread group "<<i<<" has exceeded the file size limit."<<std::endl;
		  std::cerr<<"       Rerun with more thread groups or partitions to reduce file size."<<std::endl;
		  std::cerr.flush();
		  omp_unset_lock(&token);
		  exit(2);
		}
	    }

	  ofs.close();

	  // Scrub the file name.
	  for ( j=0; j < MAX_FILE_NAME_LEN; ++j )
	    td_out_file_name[j] = '\0';

	}

      // Free up allocated memory for the next iteration, or termination of loop.
      delete [] sequences;
      delete [] bin_state;
      delete [] total_batch_length;
      
      for ( i=0; i < num_thread_groups; ++i )
	free(bins[i]);

      free(bins);
      free(bin_index);
      free(current_bin_allocation);
      free(current_batch_length);
      
    } // End OpenMP parallel for loop on the query partitions.
  
  // Free up memory at the end.
  delete [] buffer;
  delete [] query_offsets;
  delete [] q_lengths;
  delete [] qstart;

  return 0;
}
