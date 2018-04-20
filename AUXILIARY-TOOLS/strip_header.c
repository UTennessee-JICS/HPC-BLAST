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
//  Truncates the header to a reasonable length. Particularly helpful
//  for the nr/nt databases where sequence headers can be many times
//  longer than the length of short sequences.
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

void print_usage()
{
  printf("Usage:\n");
  printf("./database_distribute <options>\n");
  printf("     -h     --help     :     Prints the usage message.\n");
  printf("     -i     --input    :     REQUIRED: The name of the FASTA database to sample.\n");
  printf("     -o     --output   :     REQUIRED: The name of the output file.\n");
  printf("     -l     --length   :     REQUIRED: The length of the header to keep.\n");
  fflush(stdout);
  return;
}

int main(int argcs, char* pArgs[])
{
  // Variable declarations.

  uint64_t i = 0, y = 0, q = 0, k = 0;
  uint64_t line_length = 0;
  char * line_ptr = NULL;

  const uint64_t MAX_HEADER_LENGTH = 1024;
  const uint64_t ymax = 100;

  char ** db_seq = NULL;
  char ** seq_headers = NULL;
  uint64_t * seq_len = NULL;
  uint64_t * seq_alloc = NULL;

  uint64_t current_seq = 0;
  uint64_t num_seq_total = 0;
  uint64_t num_seq_selected = 0;
  uint64_t total_db_length = 0;
  uint64_t local_seq_length = 0;
  uint64_t line_return = 0;
  uint64_t seq_len_i=0, seq_i_bin=0;
  uint64_t seq_len_remaining = 0;
  uint64_t num_seq_alloc = 0;
  uint64_t * sample_query = NULL;
  uint64_t sample_size = 0;
  uint64_t got_newline = 0;
  uint64_t current_space = 0;

  char line[MAX_HEADER_LENGTH];

  FILE * fasta_db_in = NULL;
  FILE * fasta_db_out = NULL;

  uint64_t length=0;

  const char* const short_options = "hi:o:l:";

  const struct option long_options[] = {
    { "help", 0, NULL, 'h' },
    { "input" , 1, NULL, 'i' },
    { "output", 1, NULL, 'o' },
    { "length", 1, NULL, 'l' },
    { NULL, 0, NULL, 0 }               };

  int next_option = 0;

  int got_in = 0;
  int got_out = 0;
  int got_length= 0;

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

	case 'l':
	  length = atoi( optarg );
	  got_length = 1;
	  if ( length < 1 )
	    {
	      printf("Error: Length is non-positive.\n");
	      fflush(stdout);
	      exit(1);
	    }

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

  if ( !got_length )
    {
      length = 60;
      printf("Didn't get a length on the command line so assuming max of 60.\n");
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

  // Truncate the headers.
  for ( i=0; i < num_seq_total; ++i )
    {
      if ( strlen( seq_headers[i] ) > length )
	{
	  seq_headers[i][length-2] = '\n';
	  seq_headers[i][length-1] = '\0';
	}
    }

  for ( i=0; i < num_seq_total; ++i )
    {
      // Write the truncated sequence header.
      fprintf(fasta_db_out,"%s", seq_headers[i]);

      // Write out the amino acids.
      seq_len_remaining = seq_len[i];

      for ( y=0; y < seq_len[i]; y+=60 , seq_len_remaining-=60 )
	{
	  if ( seq_len_remaining >= 60 )
	    {
	      strncpy( line, &( db_seq[i][y] ), 60 );
	      line[60] = '\0';
	    }
	  else
	    {
	      strncpy( line, &( db_seq[i][y] ), seq_len_remaining );
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
  
  return 0;
}
