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
//  main.C
//  
//  Entry point into the program. Reads in a FASTA database file.
//  Sorts the queries by length and randomly samples them to produce
//  a query file for BLAST. Optionally, a filter on the length can be
//  applied so that only sequences of a particular lenght are sampled.
//
//  Written by - Shane Sawyer
//
//=============================================================

//  This version allows for a built-in +-1% range around the specified query length if no queries are exactly the desired length. 

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/time.h>
#include <time.h>
#include <limits.h>
#include <stdint.h>
#include <getopt.h>
#include "SFMT.h"

#define max(a,b) ({ __typeof__ (a) _a = (a); \
                    __typeof__ (b) _b = (b); \
                    _a > _b ? _a : _b; })

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

  if ( aa->length < bb->length )
    return -1;
  else if ( aa->length > bb->length )
    return 1;
  else
    return 0;
}

void print_usage()
{
  printf("Usage:\n");
  printf("./database_distribute <options>\n");
  printf("     -h     --help     :     Prints the usage message.\n");
  printf("     -i     --input    :     REQUIRED: The name of the FASTA database to sample.\n");
  printf("     -o     --output   :     REQUIRED: The name of the output file.\n");
  printf("     -s     --seqs     :     REQUIRED: The number of sequences to sample.\n");
  printf("     -f     --filter   :     Optional: A length to filter the samples by.\n");
  printf("     -b     --begin    :     Optional: The beginning of a range to filter by.\n");
  printf("     -e     --end      :     Optional: The end of a range to filter by.\n");
  printf("     -t     --total    :     Optional: The total number of letters requested by the sample.\n");

  fflush(stdout);
  return;
}

int main(int argcs, char* pArgs[])
{
  // Variable declarations.

  uint64_t i = 0, y = 0, q = 0, k = 0;
  uint64_t line_length = 0;
  uint64_t seq_min_id = 0, seq_max_id = 0;
  uint64_t is_sample_repeat = 0;

  char * line_ptr = NULL;

  const uint64_t MAX_HEADER_LENGTH = 1024;
  const uint64_t ymax = 100;
  const uint64_t SAMPLE_SIZE = 50;

  char ** db_seq = NULL;
  char ** seq_headers = NULL;
  uint64_t * seq_len = NULL;
  uint64_t * seq_alloc = NULL;
  struct node * sequences = NULL;

  //int num_buckets = 0;
  uint64_t current_seq = 0;
  uint64_t num_seq_total = 0;
  uint64_t num_seq_selected = 0;
  uint64_t total_db_length = 0;
  uint64_t local_seq_length = 0;
  uint64_t line_return = 0;
  uint64_t seq_len_i=0, seq_i_bin=0;
  uint64_t seq_min = 0, seq_avg = 0, seq_max = 0;
  double seq_avg_d = 0.0, delta_x = 0.0;
  double variance = 0.0, std_dev = 0.0;
  double s_variance = 0.0, s_std_dev = 0.0;
  double s_mean_d = 0.0;
  uint64_t s_mean = 0;
  uint64_t random_prot = 0, seed = 0;
  uint64_t seq_id = 0, seq_len_remaining = 0;
  uint64_t num_seq_alloc = 0;
  //uint64_t sample_query[SAMPLE_SIZE];
  uint64_t * sample_query = NULL;
  uint64_t sample_size = 0;
  uint64_t got_newline = 0;
  uint64_t current_space = 0;

  char line[MAX_HEADER_LENGTH];

  FILE * fasta_db_in = NULL;
  FILE * fasta_db_out = NULL;

  uint64_t filter_length=0;
  uint64_t has_filter = 0;
  uint64_t sequence_count = 0;
  uint64_t begin_range = 0;
  uint64_t end_range = 0;
  uint64_t has_range = 0;

  uint64_t target_letters = 0;

  const char* const short_options = "hi:o:s:f:b:e:t:";

  const struct option long_options[] = {
    { "help", 0, NULL, 'h' },
    { "input" , 1, NULL, 'i' },
    { "output", 1, NULL, 'o' },
    { "seqs", 1, NULL, 's' },
    { "filter", 0, NULL, 'f' },
    { "begin", 0, NULL, 'b' },
    { "end", 0, NULL, 'e' },
    { "total",0,NULL,'t'},
    { NULL, 0, NULL, 0 }               };

  int next_option = 0;

  int got_in = 0;
  int got_out = 0;
  int got_num_seqs = 0;
  int got_filter = 0;
  int got_begin = 0;
  int got_end = 0;
  int got_total = 0;

  do
    {
      next_option = getopt_long( argcs, pArgs, short_options, long_options, NULL );

      switch( next_option )
	{
	case 'h':
	  print_usage();
	  exit(0);

	case 'i':
	  fasta_db_in = fopen( optarg, "r");
	  got_in = 1;
	  
	  break;

	case 'o':
	  fasta_db_out = fopen( optarg, "w" );
	  got_out = 1;

	  break;

	case 's':
	  sample_size = atoi( optarg );
	  got_num_seqs = 1;
	  if ( sample_size < 1 )
	    {
	      printf("Error: Sample size is non-positive.\n");
	      fflush(stdout);
	      exit(1);
	    }

	  break;

	case 'f':
	  filter_length = atoi(optarg);
	  got_filter = 1;
	  has_filter = 1;
	  printf("Filtering on sequences of length %d\n", filter_length);

	  break;

	case 'b':
	  begin_range = atoi(optarg);
	  got_begin = 1;
	  has_range = 1;
	  break;

	case 'e':
	  end_range = atoi(optarg);
	  got_end = 1;
	  has_range = 1;
	  break;

	case 't':
	  target_letters = atoi(optarg);
	  got_total = 1;
	  printf("Targeting to sample up to %d total letters.\n",target_letters);

	  break;

	case '?':
	  printf("Error: Invalid command line argument!\n");
	  print_usage();
	  exit(1);

	case -1:
	  break;

	default:
	  printf("Error: Something strange happened.\n");
	  fflush(stdout);
	  exit(1);
	}

    }
  while ( next_option != -1 );
  
  if ( !got_in )
    {
      printf("Error: No input specified!\n");
      print_usage();
      exit(1);
    }

  if ( !got_out )
    {
      printf("Error: No output specified!\n");
      print_usage();
      exit(1);
    }

  // Initialize the database array to have 50,000 sequences each with an intial length of 250.
  printf("Initialize the sequence array... ");

  num_seq_alloc = 50000;
  db_seq = (char**) malloc( 50000 * sizeof(char*));

  if ( db_seq == NULL )
    {  printf("\nError: Could not allocate db_seq.\n");  exit(1);  }

  for ( i=0; i < 50000; ++i )
    {
      db_seq[i] = NULL;
      db_seq[i] = (char*) malloc( 250 * sizeof(char) );
      if ( db_seq[i] == NULL )
	{  printf("\nError: Could not allocate db_seq[%llu].\n",i);  exit(1);  }
    }

  printf("Done\n");

  // Initialize the headers for the database sequences.
  printf("Initialize the sequence header array... ");

  seq_headers = (char**) malloc( 50000 * sizeof(char*));

  if ( seq_headers == NULL )
    {  printf("\nError: Could not allocate seq_headers.\n");  exit(1);  }

  for ( i=0; i < 50000; ++i )
    {
      seq_headers[i] = NULL;
      seq_headers[i] = (char*) malloc( MAX_HEADER_LENGTH * sizeof(char) );
      if ( seq_headers[i] == NULL )
	{  printf("\nError: Could not allocate seq_headers[%llu].\n",i);  exit(1);  }
    }

  printf("Done\n");
  
  // Initialize the allocated length of the database sequences.
  printf("Initialize the sequence allocated length array... ");

  seq_alloc = (uint64_t*) malloc( 50000 * sizeof(uint64_t));

  if ( seq_alloc == NULL )
    {  printf("\nError: Could not allocate seq_alloc.\n");  exit(1);  }

  for ( i=0 ; i < 50000; ++i )
    {
      seq_alloc[i] = 250;
    }

  printf("Done\n");

  // Initialize the sequence length of the database sequences.
  printf("Initialize the sequence length array... ");

  seq_len = (uint64_t*) malloc( 50000 * sizeof(uint64_t));

  if ( seq_len == NULL )
    {  printf("\nError: Could not allocate seq_len.\n");  exit(1);  }

  for ( i=0 ; i < 50000; ++i )
    {
      seq_len[i] = 0;
    }

  printf("Done\n");

  // Read in the file one line at a time.

  printf("Reading in the database... ");

  current_seq = -1;

  if ( (line_ptr = fgets( line, MAX_HEADER_LENGTH-1, fasta_db_in )) == NULL )
    {  printf("Error: Could not get the first line of file.\n");  exit(1);  }

  while ( line_ptr != NULL )
    {
      // Check if this is a new sequence.
      if ( line[0] == '>' )  // Start of a new sequence.
	{
	  ++current_seq;

	  // Check that there is enough memory. Allocate more if needed.
	  if ( current_seq >= num_seq_alloc )
	    {
	      num_seq_alloc += (uint64_t)50000;

	      db_seq = (char**) realloc( (void*) db_seq, num_seq_alloc * sizeof(char*) );
	      if ( db_seq == NULL )
		{  printf("Error: Reallocation of db_seq failed.\n");  exit(1);  }

	      for ( i=(num_seq_alloc-50000); i < num_seq_alloc; ++i )
		{
		  db_seq[i] = NULL;
		  db_seq[i] = (char*) malloc( 250 * sizeof(char) );
		  if ( db_seq[i] == NULL )
		    {  printf("\nError: Could not allocate db_seq[%llu].\n",i);  exit(1);  }
		}

	      seq_headers = (char**) realloc( (void*) seq_headers, num_seq_alloc * sizeof(char*) );
	      if ( seq_headers == NULL )
		{  printf("Error: Reallocation of seq_headers failed.\n");  exit(1);  }

	      for ( i=(num_seq_alloc-50000); i < num_seq_alloc; ++i )
		{
		  seq_headers[i] = NULL;
		  seq_headers[i] = (char*) malloc( MAX_HEADER_LENGTH * sizeof(char) );
		  if ( seq_headers[i] == NULL )
		    {  printf("\nError: Could not allocate seq_headers[%llu].\n",i);  exit(1);  }
		}

	      seq_alloc = (uint64_t*) realloc( (void*) seq_alloc, num_seq_alloc * sizeof(uint64_t) );
	      if ( seq_alloc == NULL )
		{  printf("Error: Could not reallocate seq_alloc.\n");  exit(1);  }

	      for ( i=(num_seq_alloc-50000); i < num_seq_alloc; ++i )
		{  seq_alloc[i] = 250; }
	      
	      seq_len = (uint64_t*) realloc( (void*) seq_len, num_seq_alloc * sizeof(uint64_t) );
	      if ( seq_len == NULL )
		{  printf("Error: Could not reallocate seq_len.\n");  exit(1);  }

	      for ( i=(num_seq_alloc-50000); i < num_seq_alloc; ++i )
		{  seq_len[i] = 0; }
	    }

	  // check that we got a new line character.
	  got_newline = 0;
	  for ( i=0; i < MAX_HEADER_LENGTH; ++i )
	    {
	      if ( line[i] == '\n' )
		got_newline = 1;
	    }

	  // Copy the header information.
	  strncpy( seq_headers[current_seq], line, MAX_HEADER_LENGTH-1 );

	  current_space = MAX_HEADER_LENGTH;

	  while ( !got_newline )
	    {
	      // Grab the next line
	      line_ptr = fgets( line, MAX_HEADER_LENGTH-1, fasta_db_in );

	      // Allocate more space
	      seq_headers[current_seq] = (char*)realloc( (void*)seq_headers[current_seq], current_space+MAX_HEADER_LENGTH);

	      // Append the end.
	      strncpy( &(seq_headers[current_seq][current_space]), line, MAX_HEADER_LENGTH-1 );

	      current_space += MAX_HEADER_LENGTH;

	      for ( i=0; i < MAX_HEADER_LENGTH; ++i )
		{
		  if ( line[i] == '\n' )
		    got_newline = 1;
		}
	    }
	}
      else
	{
	  // Any other line should have nucleotides or amino acids on it. Put them into the sequence array.
	  line_length = strlen(line);

	  if ( line_length > (MAX_HEADER_LENGTH-2) )
	    {
	      printf("ERROR: Line length for sequence line is too long: %llu\n",line_length);
	      fflush(stdout);
	      exit(1);
	    }

	  // Check that the current line can be written to the sequence array. Add more memory if not.
	  if ( seq_len[current_seq] + line_length >= seq_alloc[current_seq] )
	    {
	      //printf("current_seq=%llu  seq_alloc= %llu \n",current_seq, seq_alloc[current_seq]);
	      //fflush(stdout);

	      //printf("seq_alloc[current_seq] + 250 = %llu\n", seq_alloc[current_seq] + 250 );
	      //fflush(stdout);

	      // Allocate additional space for the sequence.

	      //printf("  Realloc size requested is seq_alloc[%llu] + 250 = %llu + 250 = %llu\n",current_seq,seq_alloc[current_seq],seq_alloc[current_seq] + 250);
	      //printf("  REALLOC SIZE - (seq_alloc[current_seq] + 250 ) * (uint64_t)(sizeof(char)) = %llu \n",(seq_alloc[current_seq] + 250 ) * (uint64_t)(sizeof(char)));
	      fflush(stdout);

	      db_seq[current_seq] = (char*) realloc( (void*) db_seq[current_seq], (seq_alloc[current_seq] + 250 ) );
	      if ( db_seq[current_seq] == NULL )
		{  printf("Error: Could not reallocate db_seq[%llu].\n",current_seq);  exit(1);  }

	      seq_alloc[current_seq] += 250;
	    }

	  if ( line[line_length-1] == '\n' )
	    {  line[line_length-1] = '\0';  }

	  // Checked length already so assume its safe to pop the current line into the sequence.
	  
	  // Start from the end of the sequence ( squash the terminator already there )
	  strcpy( &(db_seq[current_seq][seq_len[current_seq]]), line );
	  
	  seq_len[current_seq] += (line_length-1);

	}

      // Read in the next line for next iteration.
      line_ptr = fgets( line, MAX_HEADER_LENGTH-1, fasta_db_in );

    }

  num_seq_total = current_seq+1;

  // Close the file.
  fclose(fasta_db_in);

  printf("Done.\n");

  printf("Total number of sequences in the database is %llu.\n", num_seq_total);

  if ( sample_size >= num_seq_total )
    {
      printf("Error: The sample size is too large. <%llu>\n",sample_size);
      fflush(stdout);
      exit(1);
    }

  sample_query = (uint64_t*) malloc(sizeof(uint64_t) * sample_size);

  if ( sample_query == NULL )
    {  printf("\nError: Could not allocate sample_query.\n");  exit(1);  }

  for ( i=0; i < sample_size; i++ )
    sample_query[i] = -1;  // Initialize to a false sequence identity.

  // Fill in the histogram array.
  
  seq_min = (uint64_t)1<<52;
  seq_max = 0;

  for ( i=0; i < num_seq_total; ++i )
    {
      seq_avg += seq_len[i];

      if ( seq_min > seq_len[i] )
	{
	  seq_min = seq_len[i];
	  seq_min_id = i;
	}

      if ( seq_max < seq_len[i] )
	{
	  seq_max = seq_len[i];
	  seq_max_id = i;
	}
    }

  seq_avg_d = ( (double)seq_avg  )/( (double) num_seq_total );

  seq_avg = (uint64_t)seq_avg_d;

  // Compute the standard deviation and variance.
  for ( i=0; i < num_seq_total; i++ )
    {
      variance += ( (double)seq_len[i] - seq_avg_d  ) * ( (double)seq_len[i] - seq_avg_d );
    }

  variance = variance * (1./( (double) (num_seq_total-1) ));

  std_dev = sqrt(variance);

  printf("Average sequence length is %lf ( or %llu ).\n", seq_avg_d, seq_avg );
  printf("Minimum sequence length is %llu.\n", seq_min);
  printf("Maximum sequence length is %llu.\n", seq_max);

  printf("\n");
  printf("Standard deviation is %lf.\n",std_dev);
  printf("Variance is %lf.\n",variance);

  if ( has_range )
    {
      if ( !got_begin )
	{
	  begin_range = seq_min;
	}
      if ( !got_end )
	{
	  end_range = seq_max;
	}
      printf("Filtering on sequences in the range of [%d , %d]\n",begin_range,end_range);
    }

  // Sort the sequences according to lengths.
  printf("Initialize the sequence node array... ");

  sequences = (struct node *) malloc( num_seq_total * sizeof(struct node));

  if ( sequences == NULL )
    {  printf("\nError: Could not allocate sequences node array.\n");  exit(1);  }

  sequence_count = 0;

  double window = 0.0;
  uint64_t lower_bound = 0;
  uint64_t upper_bound = 0;
  uint64_t midpoint = 0;
  uint64_t min_seqs = 5 * sample_size;

  if ( has_range )
    midpoint = ( end_range - begin_range )/2;

  if ( has_filter )
    {
      lower_bound = filter_length - ( (uint64_t) window * filter_length);
      upper_bound = filter_length + ( (uint64_t) window * filter_length);
    }
  else if ( has_range )
    {
      lower_bound = begin_range - ( (uint64_t) window * midpoint);
      upper_bound = end_range   + ( (uint64_t) window * midpoint);
    }
  else
    {
      lower_bound = seq_min;
      upper_bound = seq_max;
    }

  while ( sequence_count < min_seqs )
    {
      sequence_count = 0;

      for ( i=0; i < num_seq_total; ++i )
	{
	  if ( has_filter || has_range )
	    {
	      if ( seq_len[i] >= lower_bound &&
		   seq_len[i] <= upper_bound )
		{
		  sequences[sequence_count].id = i;
		  sequences[sequence_count].length = seq_len[i];
		  ++sequence_count;
		}
	    }
	  else
	    {
	      sequences[sequence_count].id = i;
	      sequences[sequence_count].length = seq_len[i];
	      ++sequence_count;
	    }
	}

      printf("Attempting to sample %llu sequences in the range of [%llu , %llu]. Found %llu in the range.\n",sample_size,lower_bound,upper_bound,sequence_count);

      if ( (has_filter || has_range) && sequence_count < min_seqs )
	{
	  window = window + 0.005;
	  lower_bound = filter_length - ( (uint64_t) (window * filter_length));
	  upper_bound = filter_length + ( (uint64_t) (window * filter_length));

	  if ( window > 0.06 ) // Ok give up.
	    {
	      fprintf(stderr,"Error, window is getting too big so I'm giving up.\n");
	      fflush(stderr);
	      exit(1);
	    }

	  printf("Adding 0.005 to the filter length window, now at %lf : ( %llu , %llu )\n",window,lower_bound,upper_bound);
	}

      if ( !(has_filter || has_range) && sequence_count == 0 )
	{
	  fprintf(stderr,"ERROR: Apparently the database is empty.\n");
	  fflush(stderr);
	  exit(1);
	}
      
    }
  printf("Done\n");

  if ( !(has_filter || has_range) )
    {
      printf("Starting the sort algorithm... ");

      // Sort the nodes by sequence length.
      qsort ( sequences , num_seq_total, sizeof(struct node), compare_length );

      printf("Done\n");
    }
  else if ( has_filter )
    {
      printf("Found %llu sequences of, or around, length %llu in the database.\n",sequence_count, filter_length);
    }
  else
    {
      printf("Found %llu sequences in the range of [%llu , %llu] in the database.\n",sequence_count,lower_bound,upper_bound);
    }

  //for (i=0;i<num_seq_total;i++)
  //  {
  //    printf("%llu     %llu\n",sequences[i].length, sequences[i].id);
  //  }

  // Initialize the psuedo-random number generator.
  //srand48(time(NULL));

  //random number = drand48();

  //mt_seed();

  //srand(time(NULL));
  seed = time(NULL);
  //init_gen_rand(1384375234u);

  //seed = 1384375254u;
  init_gen_rand(seed);

  printf("Seed used = %u\n", seed);

  uint64_t letters_collected = 0;
  uint64_t samples_taken = 0;

  for ( y=0; y < sample_size ; ++y )
    {
      is_sample_repeat = 1;
      
      while ( is_sample_repeat )
	{
	  // Assume that it is not a repeat.
	  is_sample_repeat = 0;
	  
	  // Fill the data sample.
	  random_prot = gen_rand32();
	  
	  //printf("Got random number -> %llu\n",random_prot);
	  
	  random_prot = random_prot % sequence_count;
	  
	  //printf("Modulo number of sequences -> %llu\n", random_prot);
	  
	  if ( random_prot >= sequence_count || random_prot < 0 )
	    {
	      printf("Error : random_prot = %llu\n", random_prot);
	      exit(1);
	    }
	  
	  // check if sample is a repeat.
	  for (k=0; k < y; ++k )
	    {
	      if ( random_prot == sample_query[k] )
		is_sample_repeat = 1;
	    }
	}
      
      sample_query[y] = random_prot;
      ++samples_taken;

      letters_collected += sequences[random_prot].length ;

      if ( got_total )
	{
	  if ( letters_collected > target_letters )
	    y = sample_size;  // exit early
	}
    }
  
  // possibly reset sample size
  sample_size = samples_taken;

  if ( got_total )
    {
      printf("Got a sample with %d letters and %d sequences.\n",letters_collected,sample_size);
    }

  // Compute relevant statistics.
  seq_min = (uint64_t)1<<52;
  seq_max = 0;
  seq_avg = 0;
      
  for ( y=0; y < sample_size; ++y )
    {
      q = sample_query[y];
      
      seq_avg += sequences[q].length;
      
      if ( seq_min > sequences[q].length )
	{
	  seq_min = sequences[q].length;
	  seq_min_id = q;
	}
      
      if ( seq_max < sequences[q].length )
	{
	  seq_max = sequences[q].length;
	  seq_max_id = q;
	}
    }

  seq_avg_d = ( (double)seq_avg  )/( (double) sample_size );
  
  seq_avg = (uint64_t)seq_avg_d;

  // Compute the standard deviation and variance.
  for ( y=0; y < sample_size; y++ )
    {
      q = sample_query[y];
      variance += ( (double)sequences[q].length - seq_avg_d  ) * 
	( (double)sequences[q].length - seq_avg_d );
    }
  
  variance = variance * (1./( (double) (sample_size) ));
  
  std_dev = sqrt(variance);
  
  printf("  Average sequence length is %lf ( or %llu ).\n", seq_avg_d, seq_avg );
  printf("  Minimum sequence length is %llu.\n", seq_min);
  printf("  Maximum sequence length is %llu.\n", seq_max);
  printf("  Standard deviation is %lf.\n",std_dev);
  printf("  Variance is %lf.\n\n",variance);

  for ( i=0; i < sample_size; ++i )
    {
      // Retrieve sequence id.
      seq_id = sample_query[i];
      seq_id = sequences[seq_id].id;

      // Write the sequence header.
      fprintf(fasta_db_out,"%s", seq_headers[seq_id]);

      // Write out the amino acids.
      seq_len_remaining = seq_len[seq_id];

      for ( y=0; y < seq_len[seq_id]; y+=60 , seq_len_remaining-=60 )
	{
	  if ( seq_len_remaining >= 60 )
	    {
	      strncpy( line, &( db_seq[seq_id][y] ), 60 );
	      line[60] = '\0';
	    }
	  else
	    {
	      strncpy( line, &( db_seq[seq_id][y] ), seq_len_remaining );
	      line[seq_len_remaining] = '\0';
	    }

	  fprintf(fasta_db_out,"%s\n", line);
	}

    }

  // Close the file.
  fclose(fasta_db_out);

  // Clean out the memory.
  for ( i=0; i < num_seq_total; ++i )
    {
      free(db_seq[i]);
      free(seq_headers[i]);
    }

  free(db_seq);
  free(seq_headers);
  free(seq_len);
  free(seq_alloc);
  
  free(sequences);
  free(sample_query);

  return 0;
}
