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
//  Counts up the number of letters in the query set and prints
//  out the relevant statistics.
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

#define max(a,b) ({ __typeof__ (a) _a = (a); \
                    __typeof__ (b) _b = (b); \
                    _a > _b ? _a : _b; })

struct node
{
  uint64_t id;
  uint64_t length;
};


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
  uint64_t got_newline = 0;
  uint64_t current_space = 0;

  char line[MAX_HEADER_LENGTH];

  FILE * fasta_db_in = NULL;

  if ( argcs != 2 )
    {
      printf("Incorrect commanmd line parameters specified. Just need the FASTA file.\n");
      fflush(stdout);
      exit(1);
    }

  fasta_db_in = fopen( pArgs[1], "r");

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
  printf("Total number of letters is %llu.\n",seq_avg);

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

  return 0;
}
