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
//  Entry point into the program. Reads in two output files,
//  specified on the command line, that were produced by NCBI
//  BLAST code (blast+ or blastall) and compares the results.
//  In particular, it determines the amount that the second file
//  covers from the resuls of the first file.
//
//
//  Written by - Shane Sawyer
//
//=============================================================


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/time.h>
#include <time.h>
#include <limits.h>
#include <ctype.h>
#include <getopt.h>

#define max(a,b) ({ __typeof__ (a) _a = (a); \
                    __typeof__ (b) _b = (b); \
                    _a > _b ? _a : _b; })

#define MAX_ALIGN_STR_LEN 30
#define MAX_LINE_LEN 256

struct alignment
{
  char seq_id[MAX_ALIGN_STR_LEN];
  char bit_score[MAX_ALIGN_STR_LEN];
  char evalue[MAX_ALIGN_STR_LEN];
  double score_dp;
  double evalue_dp;
};

void print_usage()
{
  printf("Usage:\n");
  printf("./compare_results <options>\n");
  printf("     -h     --help     :     Prints the usage message.\n");
  printf("     -r     --result1  :     REQUIRED: The first BLAST search file.\n");
  printf("     -R     --result2  :     REQUIRED: The second BLAST search file to compare against the first file.\n");
  printf("     -o     --out      :     Optional output file name for the comparison results.\n");

  fflush(stdout);
  return;
}

int process_blast_output( FILE* file_ptr, struct alignment ***aligns, int *num_aligns, int *num_aligns_alloced, char **query_names, int realloc_increment, char * filename );

int main(int argcs, char* pArgs[])
{
  // Variable declarations.

  int i=0, j=0, k=0;
  int num_queries=0;
  int num_queries_a=0;
  int num_queries_b=0;
  int current_query = 0;
  int current_subject = 0;

  int look_for_subjects = 0;

  int num_alignments_found_a = 0;
  int num_alignments_found_b = 0;

  char **queries_a = NULL;
  char **queries_b = NULL;

  char result_line[MAX_LINE_LEN];
  int result_line_len = 0;

  int *result_found;
  int max_results;

  char subject_id[MAX_ALIGN_STR_LEN];

  struct alignment **alignment_a = NULL;
  struct alignment **alignment_b = NULL;

  int *num_alignments_a = NULL;
  int *num_alignments_b = NULL;

  int *num_alignments_alloced_a = NULL;
  int *num_alignments_alloced_b = NULL;

  FILE *blast_file_a = NULL;
  FILE *blast_file_b = NULL;
  FILE *compare_out = NULL;

  char blast_filename_a[MAX_LINE_LEN];
  char blast_filename_b[MAX_LINE_LEN];
  char compare_outname[MAX_LINE_LEN];

  char line[MAX_LINE_LEN];

  char *line_ptr = NULL;

  char *test_line = "Query=";

  int test_line_length = 6;

  char string[MAX_LINE_LEN];
  int string_length = 0;
  int str_start_ptr=0, str_stop_ptr=0;
  int found_first_vbar = 0;
  int id_len = 0;
  int found_first_subject = 0;
  int total_alignments = 0, total_found = 0;
  int results_found = 0;

  float percent_found = 0;

  int got_blast_a = 0;
  int got_blast_b = 0;
  int got_output = 0;

  int max_seq_len = 0;
  int max_bit_len = 0;
  int max_evl_len = 0;
  int offset = 0;

  int matching_query_index = 0;

  double score_err, eval_err;
  
  const char* const short_options = "ho:r:R:";

  const struct option long_options[] = {
    { "help", 0, NULL, 'h' },
    { "out" , 1, NULL, 'o' },
    { "result1", 1, NULL, 'r' },
    { "result2", 1, NULL, 'R' },
    { NULL, 0, NULL, 0 }               };

  int next_option = 0;

  // End variable declaration.

  // Clean the empty string out.
  for ( i=0; i < MAX_LINE_LEN; i++ )
    {
      blast_filename_a[i] = '\0';
      blast_filename_b[i] = '\0';
      compare_outname[i]  = '\0';
    }

  // Parse the command line arguments.
  do
    {
      next_option = getopt_long( argcs, pArgs, short_options, long_options, NULL );

      switch( next_option )
	{
	case 'h':
	  print_usage();
	  exit(0);

	case 'o':
	  i = strlen(optarg);
	  if ( i >= MAX_LINE_LEN )
	    {
	      printf("Error: Output file name length exceeds maximum currently allowed (%d).\n",MAX_LINE_LEN);
	      exit(1);
	    }
	  strncpy(compare_outname, optarg, i);
	  got_output = 1;
	  
	  break;

	case 'r':
	  i = strlen(optarg);
	  if ( i >= MAX_LINE_LEN )
	    {
	      printf("Error: First file name length exceeds maximum currently allowed (%d).\n",MAX_LINE_LEN);
	      exit(1);
	    }
	  strncpy(blast_filename_a, optarg, i);
	  got_blast_a = 1;

	  break;

	case 'R':
	  i = strlen(optarg);
	  if ( i >= MAX_LINE_LEN )
	    {
	      printf("Error: Second file name length exceeds maximum currently allowed (%d).\n",MAX_LINE_LEN);
	      exit(1);
	    }
	  strncpy(blast_filename_b, optarg, i);
	  got_blast_b = 1;

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

  // Check that we processed the expected command line arguments.
  if ( got_blast_a == 0 )
    {
      printf("Error: The first BLAST alignment file was not specified.\n");
      print_usage();
      exit(1);
    }

  if ( got_blast_b == 0 )
    {
      printf("Error: The second BLAST alignment file was not specified.\n");
      print_usage();
      exit(1);
    }

  if ( got_output )
    {
      compare_out = fopen( compare_outname, "w");

      if ( compare_out == NULL )
	{  printf("Could not open the output file.\n");   fflush(stdout);  exit(1);  }
    }
  else
    compare_out = stdout;

  blast_file_a = fopen( blast_filename_a, "r" );
  if ( blast_file_a == NULL )
    {  printf("Could not open the first BLAST alignment file.\n");   fflush(stdout);  exit(1);  }

  blast_file_b = fopen( blast_filename_b, "r" );
  if ( blast_file_b == NULL )
    {  printf("Could not open the second BLAST alignment file.\n");   fflush(stdout);  exit(1);  }
  
  // Check that the files agree on the number of queries.

  if ( (line_ptr = fgets( line, MAX_LINE_LEN, blast_file_a )) == NULL )
    {  printf("Error: Could not get the first line of the first file.\n");  exit(1);  }

  while ( line_ptr != NULL )
    {
      if ( strncmp( test_line, line, test_line_length ) == 0 )
	{
	  // We have a query so increment the number of queries.
	  num_queries++;
	}

      // Read in the next line for next iteration.
      line_ptr = fgets( line, MAX_LINE_LEN, blast_file_a );

    }

  rewind(blast_file_a);

  if ( (line_ptr = fgets( line, MAX_LINE_LEN, blast_file_b )) == NULL )
    {  printf("Error: Could not get the first line of the second file.\n");  exit(1);  }

  while ( line_ptr != NULL )
    {
      if ( strncmp( test_line, line, test_line_length ) == 0 )
	{
	  // We have a query so decrement the number of queries.
	  num_queries_b++;
	}

      // Read in the next line for next iteration.
      line_ptr = fgets( line, MAX_LINE_LEN, blast_file_b );

    }

  rewind(blast_file_b);

  if ( (num_queries-num_queries_b) != 0 )
    {
      printf("The number of queries did not match between the two files!\n");
      printf("  The first file reported %d queries.\n", num_queries);
      printf("  The second file reported %d queries.\n", num_queries_b);
      fflush(stdout);
      fclose(blast_file_a);
      fclose(blast_file_b);
      exit(1);
    }

  // Allocate memory for the data structures.

  queries_a = (char **) malloc( sizeof(char*) * (num_queries+1) );
  if ( queries_a == NULL )
    {  printf("Could not allocate queries_a.\n");  fflush(stdout);  exit(1);  }

  for ( i=0; i <= num_queries; i++ )
    {
      queries_a[i] = (char*) malloc(sizeof(char) * MAX_ALIGN_STR_LEN );
      if ( queries_a[i] == NULL )
	{  printf("Could not allocate queries_a[%d].\n",i);  fflush(stdout);  exit(1);  }
    }

  queries_b = (char **) malloc( sizeof(char*) * (num_queries+1) );
  if ( queries_b == NULL )
    {  printf("Could not allocate queries_b.\n");  fflush(stdout);  exit(1);  }

  for ( i=0; i <= num_queries; i++ )
    {
      queries_b[i] = (char*) malloc(sizeof(char) * MAX_ALIGN_STR_LEN );
      if ( queries_b[i] == NULL )
	{  printf("Could not allocate queries_b[%d].\n",i);  fflush(stdout);  exit(1);  }
    }

  num_alignments_a = (int *) malloc( sizeof(int) * (num_queries+1) );
  if ( num_alignments_a == NULL )
    {  printf("Could not allocate num_alignments_a.\n");  fflush(stdout);  exit(1);  }

  num_alignments_b = (int *) malloc( sizeof(int) * (num_queries+1) );
  if ( num_alignments_b == NULL )
    {  printf("Could not allocate num_alignments_b.\n");  fflush(stdout);  exit(1);  }

  num_alignments_alloced_a = (int *) malloc( sizeof(int) * (num_queries+1) );
  if ( num_alignments_alloced_a == NULL )
    {  printf("Could not allocate num_alignments_alloced_a.\n");  fflush(stdout);  exit(1);  }

  num_alignments_alloced_b = (int *) malloc( sizeof(int) * (num_queries+1) );
  if ( num_alignments_alloced_b == NULL )
    {  printf("Could not allocate num_alignments_alloced_b.\n");  fflush(stdout);  exit(1);  }

  alignment_a = (struct alignment **) malloc( sizeof(struct alignment*) * ( num_queries+1 ) );
  if ( alignment_a == NULL )
    {  printf("Could not allocate alignment_a.\n");  fflush(stdout);  exit(1);  }

  alignment_b = (struct alignment **) malloc( sizeof(struct alignment*) * ( num_queries+1 ) );
  if ( alignment_b == NULL )
    {  printf("Could not allocate alignment_b.\n");  fflush(stdout);  exit(1);  }

  // Guess that each query sequence will have 10 alignments of significance.
  
  for (i=1; i <= num_queries; ++i )
    {
      num_alignments_a[i] = 0;
      num_alignments_b[i] = 0;
      
      num_alignments_alloced_a[i] = 10;
      num_alignments_alloced_b[i] = 10;

      alignment_a[i] = (struct alignment *) malloc( sizeof(struct alignment) * 10 );
      if ( alignment_a[i] == NULL )
	{  printf("Could not allocate alignment_a[%d].\n",i);  fflush(stdout);  exit(1);  }
      
      alignment_b[i] = (struct alignment *) malloc( sizeof(struct alignment) * 10 );
      if ( alignment_b[i] == NULL )
	{  printf("Could not allocate alignment_b[%d].\n",i);  fflush(stdout);  exit(1);  }

    }

  // Clean the alignment structures.
  for ( i=0; i <= num_queries; i++ )
    {
      for ( j=0; j < num_alignments_alloced_a[i]; j++ )
	{
	  for ( k=0; k < MAX_ALIGN_STR_LEN; k++ )
	    alignment_a[i][j].seq_id[k]='\0';
	  for ( k=0; k < MAX_ALIGN_STR_LEN; k++ )
	    alignment_a[i][j].bit_score[k]='\0';
	  for ( k=0; k < MAX_ALIGN_STR_LEN; k++ )
	    alignment_a[i][j].evalue[k]='\0';
	}

      for ( j=0; j < num_alignments_alloced_b[i]; j++ )
	{
	  for ( k=0; k < MAX_ALIGN_STR_LEN; k++ )
	    alignment_b[i][j].seq_id[k]='\0';
	  for ( k=0; k < MAX_ALIGN_STR_LEN; k++ )
	    alignment_b[i][j].bit_score[k]='\0';
	  for ( k=0; k < MAX_ALIGN_STR_LEN; k++ )
	    alignment_b[i][j].evalue[k]='\0';
	}
    }

  num_queries_a = process_blast_output( blast_file_a,  &alignment_a, num_alignments_a, num_alignments_alloced_a, queries_a, 10, blast_filename_a );

  fclose(blast_file_a);

  num_queries_b = process_blast_output( blast_file_b,  &alignment_b, num_alignments_b, num_alignments_alloced_b, queries_b, 10, blast_filename_b );

  fclose(blast_file_b);

  if ( num_queries_a != num_queries_b )
    {
      printf("Error: The number of queries is not consistent between the two files.\n");
      printf("  number of queries for the first file is  %d\n",num_queries_a);
      printf("  number of queries for the second file is %d\n",num_queries_b);
      fflush(stdout);
      exit(1);
    }

  // Check for running out of static memory.
  max_results = 0;
  for ( i=0; i < num_queries; i++ )
    {
      max_results = max( max_results, num_alignments_b[i] );
    }
  
  result_found = (int *) malloc( sizeof(int) * (max_results) );
  if ( result_found == NULL )
    {  printf("Could not allocate result_found. Tried to allocate %d entries.\n",max_results);  fflush(stdout);  exit(1);  }

  for ( i=1; i <= num_queries; i++ )
    {
      num_alignments_found_a += num_alignments_a[i];
      num_alignments_found_b += num_alignments_b[i];
    }

  fprintf(compare_out, "In %s, found %d alignments.\n",blast_filename_a,num_alignments_found_a);
  fprintf(compare_out, "In %s, found %d alignments.\n",blast_filename_b,num_alignments_found_b);
  fprintf(compare_out, "\n");

  // Scan the data structures looking for the longest sequence identifier, bit score, and expect value.
  max_seq_len = 4;
  max_bit_len = 12;
  max_evl_len = 8;

  for ( i=1; i <= num_queries; i++ )
    {
      for ( j=0; j < num_alignments_a[i]; j++ )
	{
	  string_length = strlen(alignment_a[i][j].seq_id);
	  max_seq_len = max( max_seq_len, string_length );

	  string_length = strlen(alignment_a[i][j].bit_score);
	  max_bit_len = max( max_bit_len, string_length );

	  string_length = strlen(alignment_a[i][j].evalue);
	  max_evl_len = max( max_evl_len, string_length );
	}
    }

  for ( i=1; i <= num_queries; i++ )
    {
      for ( j=0; j < num_alignments_b[i]; j++ )
	{
	  string_length = strlen(alignment_b[i][j].seq_id);
	  max_seq_len = max( max_seq_len, string_length );

	  string_length = strlen(alignment_b[i][j].bit_score);
	  max_bit_len = max( max_bit_len, string_length );

	  string_length = strlen(alignment_b[i][j].evalue);
	  max_evl_len = max( max_evl_len, string_length );
	}
    }

  // Compute the result line length.
  //result_line_len = max_seq_len + 2*max_bit_len + 2*max_evl_len + 18;

  result_line_len = max_seq_len + 3*max_bit_len + 3*max_evl_len + 28;

  if ( result_line_len >= MAX_LINE_LEN )
    {
      printf("ERROR: Result line is too long for static array. Need %d and have %d.\n",result_line_len,MAX_LINE_LEN);
      fflush(stdout);
      exit(1);
    }

  result_line[MAX_LINE_LEN-1] = '\0';

  // Now we can compare the results from the two files.
  for ( i=1; i <= num_queries; i++ )
    {
      // Reset the result line.
      for ( j=0; j < result_line_len; j++ )
	result_line[j]=' ';

      result_line[result_line_len] = '\0';
      
      // Reset the result_found array.
      for ( j=0; j < max_results; j++ )
	{
	  result_found[j] = 0;
	}

      // Find the index in the second file where this query occurs.
      matching_query_index = 0;
      for ( j=1; j <= num_queries; ++j )
	{
	  if ( strncmp( queries_a[i], queries_b[j], MAX_ALIGN_STR_LEN ) == 0 )
	    {
	      matching_query_index = j;
	      break;
	    }
	}

      if ( matching_query_index == 0 || j > num_queries )
	{
	  fprintf(stderr,"ERROR: For query index %d of the file input, %s, the corresponding query in the second file could not be found.\n",i,queries_a[i]);
	  fflush(stderr);
	  exit(1);
	}

      // Prep the line for blast versions used.
      string_length = strlen(blast_filename_a);
      offset = (max_bit_len+max_evl_len+7 - string_length)/2;
      strncpy(&(result_line[offset+3+max_seq_len]),blast_filename_a,string_length);

      string_length = strlen(blast_filename_b);
      offset = (max_bit_len+max_evl_len+7 - string_length)/2;
      strncpy(&(result_line[offset+8+max_seq_len+max_bit_len+max_evl_len]),blast_filename_b,string_length);

      fprintf(compare_out, "Query %d:     %s\n\n",i,queries_a[i]);

      fprintf(compare_out, "%s\n",result_line);

      // Prep the header line.
      for ( j=0; j < result_line_len; j++ )
	result_line[j]=' ';

      result_line[result_line_len] = '\0';
      result_line[max_seq_len+2] = '|';
      result_line[max_seq_len+max_bit_len+max_evl_len+12] = '|';
      result_line[max_seq_len+max_bit_len+max_evl_len+max_bit_len+max_evl_len+20] = '|';

      offset = (max_seq_len + 2 - 2)/2;
      strncpy(&(result_line[offset]),"ID",2);

      offset = (max_bit_len+4 - 9)/2;
      if ( offset < 0 ) { offset = 0; }
      strncpy(&(result_line[offset+max_seq_len+2]),"BIT SCORE",9);

      offset = (max_evl_len+3 - 6)/2;
      if ( offset < 0 ) { offset = 0; }
      strncpy(&(result_line[offset+max_seq_len+max_bit_len+3]),"EVALUE",6);

      offset = (max_bit_len+4 - 9)/2;
      if ( offset < 0 ) { offset = 0; }
      strncpy(&(result_line[offset+max_seq_len+max_bit_len+max_evl_len+10]),"BIT SCORE",9);

      offset = (max_evl_len+3 - 6)/2;
      if ( offset < 0 ) { offset = 0; }
      strncpy(&(result_line[offset+max_seq_len+2*max_bit_len+max_evl_len+13]),"EVALUE",6);

      offset = (max_bit_len+4 - 10 )/2;
      if ( offset < 0 ) { offset = 0; }
      strncpy(&(result_line[offset+max_seq_len+2*max_bit_len+2*max_evl_len+18]),"%SCORE ERR",10);

      offset = (max_evl_len+4 - 11 )/2;
      if ( offset < 0 ) { offset = 0; }
      strncpy(&(result_line[offset+max_seq_len+3*max_bit_len+2*max_evl_len+20]),"%EVALUE ERR",11);

      fprintf(compare_out,"%s\n",result_line);
      //fprintf(compare_out, "-------------------------------------------------------------\n");
      for ( j=0; j < result_line_len; j++ )
	fprintf(compare_out,"-");
      fprintf(compare_out,"\n");
      fflush(compare_out);

      // Loop over all the alignments found by the first BLAST search and print the results in that
      // order with the same output from the second BLAST search (if present).

      for ( j=0; j < num_alignments_a[i]; j++ )
	{
          // Reset the result line.
          for ( k=0; k < result_line_len; k++ )
            result_line[k]=' ';

          result_line[result_line_len] = '\0';

	  // Retrieve the subject sequence ID.
	  strcpy(subject_id, alignment_a[i][j].seq_id);

	  // Prep the first output file line output.
	  strcpy( result_line, subject_id );
	  result_line[max_seq_len+2] = '|';
	  strcpy( &(result_line[max_seq_len+4]), alignment_a[i][j].bit_score );
	  strcpy( &(result_line[max_seq_len+max_bit_len+9]), alignment_a[i][j].evalue );
	  result_line[max_seq_len+max_bit_len+max_evl_len+12] = '|';

	  for ( k=0; k < num_alignments_b[matching_query_index]; k++ )
	    {
	      if ( strcmp(subject_id , alignment_b[matching_query_index][k].seq_id) == 0 )
		{
		  result_found[k] = 1;
		  strcpy( &(result_line[max_seq_len+max_bit_len+max_evl_len+14]), alignment_b[matching_query_index][k].bit_score );
		  strcpy( &(result_line[max_seq_len+2*max_bit_len+max_evl_len+19]), alignment_b[matching_query_index][k].evalue );
		  break;
		}
	    }

	  if ( k == num_alignments_b[matching_query_index] ) // Result not found.
	    {
	      offset = (max_bit_len + max_evl_len + 6 - 16)/2;
	      if ( offset < 0 ) { offset = 1; }
	      strcpy( &(result_line[offset+max_seq_len+max_bit_len+max_evl_len+11]), "RESULT NOT FOUND" );
	    }

	  if ( k != num_alignments_b[matching_query_index] ) // Found a result.
	    {
	      // Print the relative error in bit score and evalue.
	      score_err = fabs( alignment_a[i][j].score_dp - alignment_b[matching_query_index][k].score_dp );

	      if ( fabs(alignment_a[i][j].score_dp) > 1.0e-175 )
		{
		  score_err = fabs( score_err / alignment_a[i][j].score_dp )*100.0;
		}

	      eval_err = fabs( alignment_a[i][j].evalue_dp - alignment_b[matching_query_index][k].evalue_dp );

	      if ( fabs(alignment_a[i][j].evalue_dp) > 1.0e-175 )
		{
		  eval_err = fabs( eval_err / alignment_a[i][j].evalue_dp)*100.0;
		}

	      sprintf( &(result_line[max_seq_len+2*max_bit_len+2*max_evl_len+20]), "| %.12e     %.12e",score_err, eval_err);

	    }

	  // Clean out any string terminators before the end of the string.
	  for ( k=0; k < result_line_len; k++ )
	    {
	      if ( result_line[k] == '\0' )
		result_line[k] = ' ';
	    }

	  fprintf(compare_out, "%s\n",result_line);
	}

      fprintf(compare_out, "\n");

      // Count up the number of results that were found in the second BLAST search.
      k = 0;
      results_found = 0;
      for ( j=0; j < max_results; j++ )
	{
	  if ( result_found[j] == 1 )
	    results_found++;
	}

      if ( (results_found-num_alignments_a[i]) != 0 )
	{
	  fprintf(compare_out, "Not all alignments were found in the second BLAST file.\n");
	}

      fprintf(compare_out, "Listing all other alignments in %s:\n", blast_filename_b);

      for ( j=0; j < num_alignments_b[matching_query_index]; j++ )
	{
	  if ( result_found[j] == 0 )
	    {
	      for ( k=0; k < result_line_len; k++ )
		result_line[k]=' ';
	      
	      result_line[result_line_len] = '\0';

	      strcpy( result_line, alignment_b[matching_query_index][j].seq_id );
	      result_line[max_seq_len+2] = '|';
	      strcpy( &(result_line[max_seq_len+4]), alignment_b[matching_query_index][j].bit_score );
	      strcpy( &(result_line[max_seq_len+max_bit_len+9]), alignment_b[matching_query_index][j].evalue );

	      for ( k=0; k < result_line_len; k++ )
		{
		  if ( result_line[k] == '\0' )
		    result_line[k] = ' ';
		}

	      fprintf(compare_out, "%s\n",result_line);
	    }
	}

      fprintf(compare_out, "\n");
      fprintf(compare_out, "==============================================================\n");
      
      percent_found = ( (float)results_found)/( (float)num_alignments_a[i]) * 100.0;

      fprintf(compare_out, "%s has %lf%% of the alignments for this query as %s.\n\n",blast_filename_b, percent_found, blast_filename_a);

      total_alignments += num_alignments_a[i];
      total_found += results_found;
    }

  percent_found = ( (float)total_found)/( (float)total_alignments) * 100.0;

  fprintf(compare_out, "In total, %s has %lf%% of the alignments in %s.\n",blast_filename_b, percent_found, blast_filename_a);

  
  for ( i=0; i <= num_queries; i++ )
    {
      free(alignment_a[i]);
      free(alignment_b[i]);
      free(queries_a[i]);
      free(queries_b[i]);
    }

  free(num_alignments_a);
  free(num_alignments_b);
  free(num_alignments_alloced_a);
  free(num_alignments_alloced_b);
  
  free(queries_a);
  free(queries_b);
  free(result_found);

  return 0;
}


int process_blast_output( FILE* file_ptr, struct alignment ***aligns, int *num_aligns, int *num_aligns_alloced, char **query_names, int realloc_increment, char * filename )
{
  int current_query = 0;
  int current_subject = 0;
  int i=0, j=0, k=0;
  int look_for_subjects = 0;
  int look_for_subjects_start = 0;
  char *line_ptr = NULL;
  char *test_line = "Query=";
  int test_line_length = 6;
  char line[MAX_LINE_LEN];
  char string[MAX_LINE_LEN];
  int string_length = 0;
  int str_start_ptr=0, str_stop_ptr=0;
  int found_first_vbar = 0;
  int found_second_vbar = 0;
  int id_len = 0;
  int found_first_subject = 0;
  int total_alignments = 0, total_found = 0;
  int results_found = 0;
  int num_vbar_found = 0;
  int current_char = 0;
  int found_first_char = 0;

  current_query = 0;
  current_subject = 0;
  
  if ( (line_ptr = fgets( line, MAX_LINE_LEN, file_ptr )) == NULL )
    {  printf("process_blast_output:: Error: Could not get the first line of the file. < %s >\n",filename);  exit(1);  }

  while ( line_ptr != NULL )
    {
      if ( look_for_subjects_start )
	{
	  // Two possibilities.
	  //   1:  Found alignments -- "Sequences producing...."
	  //   2:  No alignments found -- "***** No hits found *****"
	  
	  // Or keep looking for another line.
	  if ( strncmp( "Sequences producing", line, 19 ) == 0 )
	    {
	      look_for_subjects = 1;
	      look_for_subjects_start = 0;
	    }
	  else if ( strncmp( "***** No", line, 8 ) == 0 )
	    {
	      // No alignments. Move on to the next query.
	      look_for_subjects_start = 0;
	    }
	  else
	    {
	      // Do nothing. Will keep looking.
	    }

	}
      else if ( look_for_subjects )  // Flag the directs code to look for more subjects to save for the current query.
	{
	  // Scan the line and look for two vertical bars.
	  num_vbar_found = 0;

	  for ( i=0; i < MAX_LINE_LEN; ++i )
	    {
	      if ( line[i] == '|' )
		++num_vbar_found;

	      if ( line[i] == '\n' || line[i] == '\0' )
		break;
	    }

	  // At least two vertical bars indicate a subject alignment found.
	  if ( num_vbar_found >= 2 )
	    {
	      if (found_first_subject == 0)
		found_first_subject = 1;

	      // Can we accecpt another subject alignment?
	      if ( (num_aligns[current_query] + 1) >= num_aligns_alloced[current_query] )
		{
		  // Allocate more space.
		  (*aligns)[current_query] = (struct alignment*) realloc( (void*) (*aligns)[current_query],
									  ( num_aligns_alloced[current_query]+realloc_increment)
									  * sizeof(struct alignment) );

		  if ( (*aligns)[current_query] == NULL )
		    {  printf("process_blast_output:: Could not reallocate (*aligns)[%d]. < %s >\n",current_query, filename);  fflush(stdout);  exit(1);  }

		  num_aligns_alloced[current_query] = num_aligns_alloced[current_query] + realloc_increment;

		  // Clean the newly alloced space.
		  for ( j=num_aligns[current_query]; j < num_aligns_alloced[current_query]; j++ )
		    {
		      for ( k=0; k < MAX_ALIGN_STR_LEN; k++ )
			(*aligns)[current_query][j].seq_id[k]='\0';
		      for ( k=0; k < MAX_ALIGN_STR_LEN; k++ )
			(*aligns)[current_query][j].bit_score[k]='\0';
		      for ( k=0; k < MAX_ALIGN_STR_LEN; k++ )
			(*aligns)[current_query][j].evalue[k]='\0';
		    }
		}

	      // Save the aligment.

	      // First we get the subject sequence identifier, including the database identifier.

	      found_first_vbar = 0;
	      found_second_vbar = 0;
	      current_char = 0;
	      found_first_char = 0;
	      
	      while( !found_second_vbar )
		{
		  if ( line[current_char] == ' ' && !found_first_char )
		    {
		      ++current_char;
		      continue;
		    }

		  if ( !found_first_char )
		    {
		      found_first_char = 1;
		      str_start_ptr = current_char;

		      if ( line[current_char] == '|' )
			found_first_vbar = 1;

		      ++current_char;
		      continue;
		    }

		  if ( !found_first_vbar && line[current_char] == '|' )
		    {
		      found_first_vbar = 1;
		      ++current_char;
		      continue;
		    }

		  if ( found_first_vbar  && line[current_char] == '|' )
		    {
		      found_second_vbar = 1;
		      str_stop_ptr = current_char;
		    }

		  ++current_char;
		}

	      id_len = str_stop_ptr - str_start_ptr;

	      if ( id_len > MAX_ALIGN_STR_LEN )
		id_len = MAX_ALIGN_STR_LEN;

	      strncpy( &((*aligns)[current_query][current_subject].seq_id[0]) , &(line[str_start_ptr]), id_len );

	      // Now find and save the Bit Score.

	      // start the scan at position 68 of the line. find the first nonblank character.
	      // copy to the last nonblank character.
	      found_first_char = 0;

	      for ( i=68; i < MAX_LINE_LEN; ++i )
		{
		  if ( !found_first_char )
		    {
		      if ( line[i] == ' ' )
			continue;
		      else
			{
			  found_first_char = 1;
			  str_start_ptr = i;
			}
		    }

		  if ( line[i] == ' ' )
		    {
		      break;
		    }

		  if ( line[i] == '\0' )
		    break;
		}

	      if ( i >= 150 || (i-str_start_ptr) >= MAX_ALIGN_STR_LEN )  // didn't find the Bit Score.
		{
		  printf("process_blast_output:: Encountered an error in getting the Bit Score from the file. < %s >\n", filename);
		  printf("process_blast_output:: position of end of Bit Score = %d,  query number = %d\n",i,current_query);
		  fflush(stdout);
		  exit(1);
		}

	      strncpy( &((*aligns)[current_query][current_subject].bit_score[0]) , &(line[str_start_ptr]), (i-str_start_ptr) );

	      (*aligns)[current_query][current_subject].score_dp = atof( &((*aligns)[current_query][current_subject].bit_score[0]) );
	      
	      // Now find and save the Evalue.

	      // start looking from where the bit score left off.

	      found_first_char = 0;

	      for ( ; i < MAX_LINE_LEN; ++i )
		{
		  if ( !found_first_char )
		    {
		      if ( line[i] == ' ' )
			continue;
		      else
			{
			  found_first_char = 1;
			  str_start_ptr = i;
			}
		    }

		  if ( line[i] == ' ' )
		    {
		      break;
		    }

		  if ( line[i] == '\0' || line[i] == '\n' )
		    {
		      break;
		    }
		}

	      if ( i >= 180 || (i-str_start_ptr) >= MAX_ALIGN_STR_LEN )  // didn't find the Evalue.
		{
		  printf("process_blast_output:: Encountered an error in getting the Expect Value from the file. < %s >\n", filename);
		  printf("process_blast_output:: position of end of Expect Value = %d,  query number = %d\n",i,current_query);
		  fflush(stdout);
		  exit(1);
		}

	      strncpy( &((*aligns)[current_query][current_subject].evalue[0]) , &(line[str_start_ptr]), (i-str_start_ptr) );

	      (*aligns)[current_query][current_subject].evalue_dp = atof( &((*aligns)[current_query][current_subject].evalue[0]) );

	      num_aligns[current_query] = num_aligns[current_query] + 1;

	      current_subject++;

	    }
	  else if ( !found_first_subject )
	    {

	    }
	  else  // No more subjects.
	    {
	      look_for_subjects = 0;
	    }
	}
      else if ( strncmp( test_line, line, test_line_length ) == 0 )
	{
	  look_for_subjects_start = 1;

	  found_first_subject = 0;

	  current_subject = 0;

	  current_query++;

	  if ( query_names != NULL )  	      // Save the query id.
	    {
	      // Query name will include the database identifier and the identifier within the vertical bars.

	      //Start at just after 'Query='
	      found_first_vbar = 0;
	      found_second_vbar = 0;
	      current_char = 6;
	      found_first_char = 0;

	      while( !found_second_vbar )
		{
		  if ( line[current_char] == ' ' && !found_first_char )
		    {
		      ++current_char;
		      continue;
		    }

		  if ( !found_first_char )
		    {
		      found_first_char = 1;
		      str_start_ptr = current_char;

		      if ( line[current_char] == '|' )
			found_first_vbar = 1;

		      ++current_char;
		      continue;
		    }

		  if ( !found_first_vbar && line[current_char] == '|' )
		    {
		      found_first_vbar = 1;
		      ++current_char;
		      continue;
		    }

		  if ( found_first_vbar  && line[current_char] == '|' )
		    {
		      found_second_vbar = 1;
		      str_stop_ptr = current_char;
		    }

		  ++current_char;
		}

	      id_len = str_stop_ptr - str_start_ptr;

	      if ( id_len > MAX_ALIGN_STR_LEN )
		id_len = MAX_ALIGN_STR_LEN;

	      strncpy(query_names[current_query],  &(line[str_start_ptr]), id_len );
	    }

	}

      // Read in the next line for next iteration.
      line_ptr = fgets( line, MAX_LINE_LEN, file_ptr );

    }

  return current_query;
}
