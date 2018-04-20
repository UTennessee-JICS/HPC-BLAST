/*  Copyright 2016 UTK JICS AACE
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *           http://www.apache.org/licenses/LICENSE-2.0
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 *
 *
 * ===========================================================================
 *  $Id: blastx_app.cpp 461340 2015-03-09 18:08:15Z ivanov $
 * ===========================================================================
 *
 *                            PUBLIC DOMAIN NOTICE
 *               National Center for Biotechnology Information
 *
 *  This software/database is a "United States Government Work" under the
 *  terms of the United States Copyright Act.  It was written as part of
 *  the author's official duties as a United States Government employee and
 *  thus cannot be copyrighted.  This software/database is freely available
 *  to the public for use. The National Library of Medicine and the U.S.
 *  Government have not placed any restriction on its use or reproduction.
 *
 *  Although all reasonable efforts have been taken to ensure the accuracy
 *  and reliability of the software and data, the NLM and the U.S.
 *  Government do not and cannot warrant the performance or results that
 *  may be obtained by using this software or data. The NLM and the U.S.
 *  Government disclaim all warranties, express or implied, including
 *  warranties of performance, merchantability or fitness for any particular
 *  purpose.
 *
 *  Please cite the author in any work or product based on this material.
 *
 * ===========================================================================
 *
 * Authors:  Christiam Camacho
 *
 */

/** @file hpc_blastx_app.cpp
 * BLASTX command line application
 */

#ifndef SKIP_DOXYGEN_PROCESSING
static char const rcsid[] = 
	"$Id: blastp_app.cpp 461340 2015-03-09 18:08:15Z ivanov $";
#endif /* SKIP_DOXYGEN_PROCESSING */

#include <mpi.h>

#include <string>
#include <iostream>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <stdint.h>

#include <ncbi_pch.hpp>
#include <corelib/ncbiapp.hpp>
#include <algo/blast/api/local_blast.hpp>
#include <algo/blast/api/remote_blast.hpp>
#include <algo/blast/blastinput/blast_fasta_input.hpp>
#include <algo/blast/blastinput/blastx_args.hpp>
#include <algo/blast/api/objmgr_query_data.hpp>
#include <algo/blast/format/blast_format.hpp>
#include "blast_app_util.hpp"

#include "../../algo/blast/api/blast_aux_priv.hpp"
#include "../../algo/blast/api/blast_seqalign.hpp"
#include "../../algo/blast/core/blast_psi_priv.h"//for serialization
#include <objects/seq/Seq_annot.hpp>//for output of alignments ceb
#include <omp.h>
#include <fstream>

#ifndef SKIP_DOXYGEN_PROCESSING
USING_NCBI_SCOPE;
USING_SCOPE(blast);
USING_SCOPE(objects);
#endif

#define USE_RESTART 0
#define UPDATE 1
#define TIMING 0
#define OUTPUT 1
#define CLEAN_UP_OBJECTS 1
#define PARTIAL_DATABASE 0

// HPC-BLAST::
// Define some global variables for MPI.

int g_rank, g_size;
int rep_rank, rep_group;

// Define globals for thread parameters
int num_thread_groups = 0;
int num_team_leaders = 0;
int num_searcher_threads = 0;
int num_threads_total = 0;
int num_ranks_per_group = 0;
int num_replication_groups = 0;
int num_ranks_total = 0;
std::string num_threads_str;
std::string db_str;
std::string query_str;
std::string out_str;


inline bool file_exists(const std::string& name)
{
  struct stat buffer;   
  return (stat (name.c_str(), &buffer) == 0); 
}

struct membuf : std::streambuf
{
  membuf(char *buf, char *end) { this->setg(buf,buf,end); }
};

// Combine HSPResults within a thread group (across team leaders)
// // Input is array of serialized HSPResults objects
BlastHSPResults* gbhrpfncCombineTeamLeaderBlastHSPResults  (char ** HSPResultsBuffers,
                                                            size_t * buffer_size,
                                                            int tid,
                                                            int num_team_leaders );

// Combine BlastHSPResults across team leaders
// Input is array of HSPResults objects
BlastHSPResults* gbhrpfncCombineTeamLeaderBlastHSPResults ( BlastHSPResults** abhrppBlastHSPResults,
                                                            int aintThreadID,
                                                            int aintNumberOfTeamLeaders );

// Output objects to file
void gvoifncDumpBlastHSPResults (BlastHSPResults *abhrpBlastHSPResults);

void gvoifncDumpBlastQueryInfo (BlastQueryInfo *abqipBlastQueryInfo);

void gvoifncDumpBlast_KarlinBlk (Blast_KarlinBlk **abkbppBlast_KarlinBlk,
                                 int aintNumberOfContexts);

void gvoifncDumpBlast_KarlinBlk1 (Blast_KarlinBlk *abkbpBlast_KarlinBlk);

void gvoifncDumpSBlastScoreMatrix(SBlastScoreMatrix *absmpSBlastScoreMatrix,int alphabet_size);

void gvoifncDumpBlastScoreBlk(BlastScoreBlk *absbpBlastScoreBlk);

//Combine objects across thread groups
BlastScoreBlk *gbsbpfncCombineThreadGroupBlastScoreBlk(int aintNumberOfThreadGroups,
                                                       BlastScoreBlk **ScoreBlkGroup);

BlastHSPResults *gbhrpfncCombineThreadGroupBlastHSPResults(int aintNumberOfThreadGroups,
                                                           BlastHSPResults **abhrppBlastHSPResults);

BlastQueryInfo *gbqipfncCombineThreadGroupBlastQueryInfo(int aintNumberOfThreadGroups,
                                                         BlastQueryInfo **abqippBlastQueryInfo);

CRef<CBlastQueryVector> gbqvfncCombineReplicationGroupQueries(int aintNumberOfThreadGroups,
                                                              CRef<CBlastQueryVector>* abqvpQueryBatchGroup);

BlastHSPResults* gbhrpfncDuplicateBlastHSPResults( BlastHSPResults* HSPResultsArray);
BlastScoreBlk *gbsbpfncDuplicateBlastScoreBlk(BlastScoreBlk *ScoreBlkGroup);


CRef<CBlastQueryVector> gbqvfncCombineReplicationGroupQueries(int aintNumberOfThreadGroups,
                                                              char** buffer,
                                                              int gbl_bsize,
                                                              int* lintpBatchStartPosition,
                                                              int* lintpBatchEndPosition,
                                                              CBlastInputSourceConfig* iconfig,
                                                              int QueryBatchSize,
                                                              int gbl_query_count,
                                                              CRef<CScope> gbl_scope);

//Create search results set from objects
CRef<CSearchResultSet> gsrsfncCollectResults(BlastHSPResults* hsp_results,
                                             BlastQueryInfo* m_QueryInfo,
                                             BlastScoreBlk* m_ScoreBlk,
                                             CRef<IQueryFactory> m_QueryFactory,
                                             const CBlastOptions* m_Options,
                                             CRef<IBlastSeqInfoSrc> m_SeqInfoSrc);


BlastHSPList* BlastHSPListDuplicate(const BlastHSPList *bhspl);

//==========================================




class CBlastxApp : public CNcbiApplication
{
public:
    /** @inheritDoc */
    CBlastxApp() {
        CRef<CVersion> version(new CVersion());
        version->SetVersionInfo(new CBlastVersion());
        SetFullVersion(version);
    }
private:
    /** @inheritDoc */
    virtual void Init();
    /** @inheritDoc */
    virtual int Run();

    /// This application's command line args
    CRef<CBlastxAppArgs> m_CmdLineArgs;

   //ceb array of pointers for each threads objects necessary to build output 
   BlastHSPResults ** tHSPResults = NULL;
   BlastQueryInfo ** tQueryInfo = NULL;
   BlastScoreBlk ** tBlastScoreBlk = NULL;

   BlastHSPResults ** HSPResultsGroup = NULL;
   BlastQueryInfo ** QueryInfoGroup = NULL;
   BlastScoreBlk ** BlastScoreBlkGroup = NULL;
   CRef<CBlastQueryVector>* QueryBatchGroup = NULL;
};

//ceb function to test for existence of file
bool fexists(const char *filename)
{
  ifstream ifile(filename);
  return ifile.good();
};


//ceb This function searches a file for n occurances of a given phrase.
//    It returns the file position of the line following the 
//    nth phrase location. If n occurances of the phrase are not found,
//    the NULL value is returned.
long get_last_good_position(string filename, string phrase, int nqueries)
{
    std::ifstream infile(filename,std::fstream::in);
    if ( !infile.is_open() ) {
        cerr << "Error while opening file." << endl;
        exit( EXIT_FAILURE );
    }    

    int count=0;
    string line;   
    //printf("nqueries=%d\n",nqueries);
    //printf("searching output file for phrase: %s\n",phrase.c_str());
    while( !infile.eof() ) 
    {
       std::getline(infile,line);  
       //printf("current line: %s\n",line.c_str());
       if( strstr(line.c_str(), phrase.c_str()) != NULL)
       {
           count++;
           //printf("%dth phrase %s found\n",count,phrase.c_str());
           if (count==nqueries) 
           {
               long pos = infile.tellg();
               //printf("last completed report ended at pos %d\n",pos);
               infile.close();
               //return current position in the file
               return pos;
           }
       }
    }   
    //printf("%dth phrase %s found\n",count,phrase.c_str());
    infile.close();
    return NULL;
}


void CBlastxApp::Init()
{
    // formulate command line arguments

    m_CmdLineArgs.Reset(new CBlastxAppArgs());

    // read the command line

    HideStdArgs(fHideLogfile | fHideConffile | fHideFullVersion | fHideXmlHelp | fHideDryRun);
    SetupArgDescriptions(m_CmdLineArgs->SetCommandLine());
}

int CBlastxApp::Run(void)
{
#if TIMING
    double t_timer[100]={0.};
    t_timer[1] = omp_get_wtime();
#endif
    int status = BLAST_EXIT_SUCCESS;

    bool working = true;

    // Declare the rank's buffer.
    char **buffer = 0;
    int *bsize = 0;

    int n_batches_precomputed=0;


    buffer = (char**) malloc( num_thread_groups * sizeof(char*) );
    if ( buffer == 0 ){ return -1; }
    for ( int i=0; i < num_thread_groups; ++i ){ buffer[i] = 0; }

    bsize = (int*) malloc( num_thread_groups * sizeof(int) );
    if ( bsize == 0 ){ return -1; }

    // Query buffer end posititon for each query by thread group
    int **lintppThreadGroupQueryBufferEndPosition=0;
    lintppThreadGroupQueryBufferEndPosition=(int**)malloc(num_thread_groups*sizeof(int*));
    if (lintppThreadGroupQueryBufferEndPosition==0) {return -1;}
    
    // Query buffer batch start position by thread group
    int *lintpBatchStartPosition=0;
    lintpBatchStartPosition=(int*)malloc(num_thread_groups*sizeof(int));
    if (lintpBatchStartPosition==0) {return -1;}
    
    // Query buffer batch end position by thread group
    int *lintpBatchEndPosition=0;
    lintpBatchEndPosition=(int*)malloc(num_thread_groups*sizeof(int));
    if (lintpBatchEndPosition==0) {return -1;}


    // Write out the global vars.
    std::cout<<"COMM_WORLD has "<<g_size<<" ranks."<<std::endl;
    std::cout<<"My rank is "<<g_rank<<std::endl;
    std::cout<<"My rep_group is "<<rep_group<<" and rep_rank is "<<rep_rank<<std::endl;
    std::cout<<"num_thread_groups = "<<num_thread_groups<<std::endl;
    std::cout<<"num_team_leaders = "<<num_team_leaders<<std::endl;
    std::cout<<"num_searcher_threads = "<<num_searcher_threads<<std::endl;
    std::cout<<"num_ranks_per_group = "<<num_ranks_per_group<<std::endl;
    std::cout<<"num_replication_groups = "<<num_replication_groups<<std::endl;

    int gbl_bsize=0;//ceb

    //Load query file into buffer
    for ( int i=0; i < num_thread_groups; ++i )
      {

	// Get the file name.
	std::stringstream qfile_name;
	qfile_name << query_str << "." <<rep_group<< "." <<i;

	// Open a stream to the file.
	std::ifstream ifs(qfile_name.str().c_str(), std::ifstream::binary );

	// Get the size of the file.
	ifs.seekg(0,ifs.end);
	bsize[i] = ifs.tellg();
	ifs.seekg(0,ifs.beg);  // Go back to the beginning.
	
	// Allocate the buffer.
	buffer[i] = (char*) malloc( (bsize[i]+1) * sizeof(char) );

        if ( buffer[i] == 0 ){ return -1; }

	// Read the file into the buffer.
	ifs.read( buffer[i], bsize[i] );

	// Close the stream.
	ifs.close();
#if 0
	std::cout<<"Read in "<<bsize[i]<<" bytes from file "<<qfile_name.str()<<std::endl;
#endif

       gbl_bsize+=bsize[i];


       //Set up variables for use in tracking postion of input buffers for output merging 
       // Find number of queries in this thread group
       int lintNumberOfThreadGroupQueries=0;
       for (int j=0; j<bsize[i]; ++j) {
	 //a query entry is denoted by the ">gi" sequence of characters
	 if (buffer[i][j]=='>') {  
           //a query entry can have multiple >gi entries prior to the listing of the sequence
	   // we need to wait till we reach a return character before we can count the query
	   while(buffer[i][j]!='\n' && j<bsize[i]){
	     ++j;
	   }
	   lintNumberOfThreadGroupQueries++; 
	 }
       }


       //std::cout<<"lintNumberOfThreadGroupQueries="
       //<<lintNumberOfThreadGroupQueries<<std::endl;std::cout.flush();

       // Allocate space for query buffer end position by thread group
       lintppThreadGroupQueryBufferEndPosition[i]=(int*)malloc(lintNumberOfThreadGroupQueries*sizeof(int));
       
       // Set query buffer end position by thread group
       int lintCurrentQuery=-1;
       for (int j=0; j<bsize[i]; ++j) {
	 //a query entry is denoted by the ">gi" sequence of characters
	 if (buffer[i][j]=='>'){
	   //the next query must be preceded by one or more return characters
	   
	   while(j<bsize[i])
	   {
	     ++j;
	     if(buffer[i][j]=='\n' && buffer[i][j+1]=='>'){break;}
	   }

	   lintCurrentQuery++;
	   if (j>1) {
	     lintppThreadGroupQueryBufferEndPosition[i][lintCurrentQuery]=j;
	     //std::cout<<"lintppThreadGroupQueryBufferEndPosition["<<i<<"]["<<lintCurrentQuery<<"]="<<lintppThreadGroupQueryBufferEndPosition[i][lintCurrentQuery-1]<<std::endl;std::cout.flush();
	   }
	 }

       }

       // Set final query buffer end position 
       lintppThreadGroupQueryBufferEndPosition[i][lintNumberOfThreadGroupQueries-1]=bsize[i]-1;
       //std::cout<<"lintppThreadGroupQueryBufferEndPosition["<<i<<"]["
       //<<lintNumberOfThreadGroupQueries-1<<"]="
       //<<lintppThreadGroupQueryBufferEndPosition[i][lintNumberOfThreadGroupQueries-1]
       //<<std::endl;std::cout.flush();
      }
    std::cout.flush();

    // Allocate the thread buffers for serialized data (only the first ptr, the threads will allocate for us)
    int number_of_threads = num_thread_groups*num_team_leaders;

    // Pointers to thread objects used to consolidate output
    tHSPResults = new BlastHSPResults*[number_of_threads];
    tQueryInfo = new BlastQueryInfo*[number_of_threads];
    tBlastScoreBlk = new BlastScoreBlk*[number_of_threads];

    //variables to point to objects merged within each thread group
    HSPResultsGroup = new BlastHSPResults*[num_thread_groups];
    QueryInfoGroup = new BlastQueryInfo*[num_thread_groups];
    BlastScoreBlkGroup = new BlastScoreBlk*[num_thread_groups];
    QueryBatchGroup = new CRef<CBlastQueryVector>[num_thread_groups];

    std::cout << "Address of buffer is " << buffer << std::endl;
    for ( int tg=0; tg < num_thread_groups; tg++ )
      std::cout << "Address of buffer["<< tg << "] is " << &(buffer[tg]) << std::endl;
    std::cout << "Address of bsize is " << bsize << std::endl;
    for ( int tg=0; tg < num_thread_groups; tg++ )  {
      std::cout << "Address of bsize[" << tg << "] is " << &(bsize[tg]) << std::endl;
      std::cout << "Contents of bsize[" << tg << "] is " << bsize[tg] << std::endl;
    }
    std::cout << "Value of gbl_bsize is " << gbl_bsize << std::endl;
    std::cout.flush();

#pragma omp parallel default(shared) num_threads( num_thread_groups*num_team_leaders )
    {
      int tid = omp_get_thread_num();

      // Compute the thread group I belong to.
      int thread_group_id = tid / num_team_leaders;
      int local_tid = tid % num_team_leaders;

       // Running total of thread group queries emitted from GetNextSeqBatch
       int lintCurrentNumberOfThreadGroupQueries=0;
       
       // Calculate the next query buffer start position from the current query buffer end position (add 1) 
       int lintNextBatchStartPosition=0;

       //ceb
       int gbl_query_count=0;

      try 
      {
        // Allow the fasta reader to complain on invalid sequence input
        SetDiagPostLevel(eDiag_Warning);
        SetDiagPostPrefix("hpc_blastx");

        /*** Get the BLAST options ***/
        const CArgs& args = GetArgs();

	CArgs my_args(GetArgs());

	// try to initialize a local copy.
	CRef<CBlastxAppArgs> my_CmdLineArgs;
	my_CmdLineArgs.Reset(new CBlastxAppArgs());

        CRef<CBlastOptionsHandle> opts_hndl;
        if(RecoverSearchStrategy(my_args, my_CmdLineArgs)) 
        {
        	opts_hndl.Reset(&*my_CmdLineArgs->SetOptionsForSavedStrategy(my_args));
        }
        else 
        {
           opts_hndl.Reset(&*my_CmdLineArgs->SetOptions(my_args));
        }

        const CBlastOptions& opt = opts_hndl->GetOptions();

        /*** Initialize the database/subject ***/
        CRef<CBlastDatabaseArgs> db_args(my_CmdLineArgs->GetBlastDatabaseArgs());

	CRef<CSearchDatabase> my_search_db = db_args->GetSearchDatabase();
	CRef<CSeqDB> my_seq_db = my_search_db->GetSeqDb();

	// Retreive the extents of the database from the CSeqDB object and compute the
	// actual extents threads in a given thread group will use.
	int db_num_seqs = my_seq_db->GetNumSeqs();

	int num_seqs_per_team_leader, my_oid_start, my_oid_stop;


	if ( db_num_seqs < num_team_leaders )
	  {
	    // Assign each available team leader at most a single database sequence.
	    if ( local_tid < db_num_seqs )
	      {
		my_oid_start = (local_tid + 1);
		my_oid_stop =  (local_tid + 2);
	      }
	    else
	      {
		my_oid_start = 0;
		my_oid_stop  = 0;
	      }
	  }
	else
	  {
	    num_seqs_per_team_leader = (db_num_seqs + num_team_leaders-1)/num_team_leaders;
	    my_oid_start =  local_tid      *num_seqs_per_team_leader;
	    my_oid_stop = ( local_tid + 1 )*num_seqs_per_team_leader;
	     if ( local_tid == (num_team_leaders-1) ){  my_oid_stop = db_num_seqs; }
	  }

//ceb
#if PARTIAL_DATABASE
my_oid_stop =  ( local_tid + 1 )*(num_seqs_per_team_leader/10);
#endif

	// Reset this thread's view into the database.
	my_seq_db->SetIterationRange(my_oid_start, my_oid_stop);

        CRef<CLocalDbAdapter> db_adapter;
#if UPDATE
        CRef<CScope> lcl_scope;//ceb
        InitializeSubject(db_args, opts_hndl, my_CmdLineArgs->ExecuteRemotely(),
			  db_adapter, lcl_scope);//ceb
        _ASSERT(db_adapter && lcl_scope);

	//create a global scope object for thread 0
	CRef<CScope> gbl_scope;
	if(tid==0) InitializeSubject(db_args, opts_hndl, my_CmdLineArgs->ExecuteRemotely(),db_adapter, gbl_scope);//ceb
#else
        CRef<CScope> scope;
        InitializeSubject(db_args, opts_hndl, my_CmdLineArgs->ExecuteRemotely(),
                          db_adapter, scope);
        _ASSERT(db_adapter && scope);
#endif

	#pragma omp barrier

        /*** Get the query sequence(s) ***/
        CRef<CQueryOptionsArgs> query_opts = my_CmdLineArgs->GetQueryOptionsArgs();

        SDataLoaderConfig dlconfig =
            InitializeQueryDataLoaderConfiguration(query_opts->QueryIsProtein(),
                                                   db_adapter);

        CBlastInputSourceConfig iconfig(dlconfig, query_opts->GetStrand(),
                                     query_opts->UseLowercaseMasks(),
                                     query_opts->GetParseDeflines(),
                                     query_opts->GetRange());
#if UPDATE
	//ceb
	  //Set up thread local input stream
	  //buffer[i] contains full query input from each file[i]
	  membuf lcl_query_membuf( &(buffer[thread_group_id][0]),
				   &(buffer[thread_group_id][ bsize[thread_group_id] - 1 ]) );
	  auto_ptr<CNcbiIstream> lcl_query_input_stream;
	  lcl_query_input_stream.reset(new CNcbiIstream(&lcl_query_membuf));
          if(IsIStreamEmpty(*lcl_query_input_stream)){ERR_POST(Warning << "Query is Empty!");}
          CBlastFastaInputSource lcl_fasta(*lcl_query_input_stream, iconfig);
          CBlastInput lcl_input(&lcl_fasta, my_CmdLineArgs->GetQueryBatchSize());
#else
	membuf my_membuf( &(buffer[thread_group_id][0]), 
	  	  &(buffer[thread_group_id][ bsize[thread_group_id] - 1 ]) );
	auto_ptr<CNcbiIstream> my_new_input_stream;
	my_new_input_stream.reset(new CNcbiIstream(&my_membuf));
        if(IsIStreamEmpty(*my_new_input_stream)){
             	ERR_POST(Warning << "Query is Empty!");
	//return BLAST_EXIT_SUCCESS;
	}
        CBlastFastaInputSource fasta(*my_new_input_stream, iconfig);
        CBlastInput input(&fasta, my_CmdLineArgs->GetQueryBatchSize());
#endif

	// Redirect thread's output file.
	string my_string (args[kArgOutput].AsString());
	string my_output;
	std::stringstream sstream;

	// Construct the thread's output file name.
	sstream << my_string << "." << rep_group << "." << rep_rank << "." << thread_group_id << "." << local_tid;
	my_output = sstream.str();
	
	auto_ptr<CNcbiOstream> my_new_output;
#if UPDATE
#if USE_RESTART //Restart code

	  //ceb Check to see if restart file exists for this output file (set of queries)
	  //    If restart file does not exist, then we default to no restart.
	  //    If we fail gracefully at the end, we will create a restart file
	  //    Note: restart files correspond to thread output, so there may be restart files for some threads
	  //    and not for others.
	  bool restarted=false;
	  int my_batches_precomputed=0;
	  int query_count=0;
	  string my_restart;
	  auto_ptr<CNcbiFstream> restart_stream;
	  if(tid==0)
	    {
	      //ceb create restart file name. Identical to output file name but with .restart appended
	      sstream << ".restart";
	      my_restart = sstream.str();

	      //ceb Open stream for restart file
	      restart_stream.reset(new CNcbiFstream(my_restart.c_str(),std::ios::in));
	      if(restart_stream->good())
		{  //ceb If file exists, then open input stream and read count
		  *restart_stream >> my_batches_precomputed >> query_count;
		  //printf("Opening restart file for reading: %d batches with %d queries precomputed\n",n_batches_precomputed,query_count);
		  restart_stream->close();//close this file. We may recreate it at the end
		  if(my_batches_precomputed>0) restarted=true;
		}


              //ceb Setting up output file
              if(restarted && fexists(my_output.c_str()))
		{ //ceb If file exists, then open output stream and begin write at end of file
		  //printf("Opening previous output file\n");
		  //ceb locate pointer location of last good position based on number of batch queries successfully written
		  string phrase = "Effective search space used:";//phrase indicates end of query output  
		  long pos = get_last_good_position(my_output.c_str(), phrase, query_count);
		  //printf("restart file pos=%d\n",pos);
		  //ceb open file as output and move write pointer to position
		  //File must be read/write in order to alter the position of the initial write
		  CNcbiOfstream* my_output_stream = new CNcbiOfstream(my_output.c_str(),std::fstream::in | std::fstream::out);
		  //move pointer
		  if(pos!=NULL) 
		    {
		      my_output_stream->seekp(pos,ios_base::beg);
		    }
		  //assign pointer to stream
		  my_new_output.reset(my_output_stream);
		}
	      else
		{  //ceb else create new file and open for writing
		  //printf("Opening new output file\n");
		  my_new_output.reset(new CNcbiOfstream(my_output.c_str()));
		  //ceb if a restart file exists but the output file does not, we have to start from scratch
		  restarted=false;
		  my_batches_precomputed=0;
		}
	    } 


	  //make sure all threads get the number of precomputed batches
          //#pragma omp single
	    //if(n_batches_precomputed < my_batches_precomputed)
	  //{n_batches_precomputed = my_batches_precomputed;}
          if (tid==0){n_batches_precomputed = my_batches_precomputed;}
          #pragma omp flush(n_batches_precomputed)
          #pragma omp barrier

#else //no restart
	  if (tid==0){ my_new_output.reset(new CNcbiOfstream(my_output.c_str()));}
#endif//restart
#else
	  {my_new_output.reset(new CNcbiOfstream(my_output.c_str()));}
#endif//UPDATE

        /*** Get the formatting options ***/
        CRef<CFormattingArgs> fmt_args(my_CmdLineArgs->GetFormattingArgs());
	CBlastFormat* formatter;//ceb
#if UPDATE
        if(tid==0)
#endif
	{
        formatter = new CBlastFormat(opt, *db_adapter,
                               fmt_args->GetFormattedOutputChoice(),
                               query_opts->GetParseDeflines(),
                               *my_new_output,
                               fmt_args->GetNumDescriptions(),
                               fmt_args->GetNumAlignments(),
#if UPDATE
				     *gbl_scope,
#else
				     *scope,
#endif
                               opt.GetMatrixName(),
                               fmt_args->ShowGis(),
                               fmt_args->DisplayHtmlOutput(),
                               opt.GetQueryGeneticCode(),
                               opt.GetDbGeneticCode(),
                               opt.GetSumStatisticsMode(),
                               my_CmdLineArgs->ExecuteRemotely(),
                               db_adapter->GetFilteringAlgorithm(),
                               fmt_args->GetCustomOutputFormatSpec());
        
        formatter->SetQueryRange(query_opts->GetRange());
        //formatter->SetLineLength(fmt_args->GetLineLength());
        if((fmt_args->GetFormattedOutputChoice() ==  CFormattingArgs::eXml2 ||
            fmt_args->GetFormattedOutputChoice() ==  CFormattingArgs::eJson)
           && args[kArgOutput].AsString() != "-")
	  {formatter->SetBaseFile(args[kArgOutput].AsString());}
#if USE_RESTART
            //ceb Skip printing the file prologue if we are restarting from a previous solution
            if(!restarted){formatter->PrintProlog();}
#else
//#if OUTPUT
        formatter->PrintProlog();
//#endif
#endif //RESTART
	}
        #pragma omp single //default thread global working var to TRUE
	{working = true;}

	bool last_search = false;
	bool got_work = false;

	#pragma omp barrier //wait for global working variable to be set

        /*** Process the input ***/
	int batch_count=0;

#if UPDATE
        for(; working; )
#else
        for (; working; formatter->ResetScopeHistory())
#endif
         {
#if UPDATE
	    //clear memory from last scope object
	    if(tid==0){formatter->ResetScopeHistory();}
	    lcl_scope->ResetDataAndHistory();//ceb
#endif //UPDATE

	  my_seq_db->SetIterationRange(my_oid_start, my_oid_stop);
	  #pragma omp barrier //wait for window into db to be set (same for all threads in rank)
#if UPDATE
          CRef<CBlastQueryVector> query_batch(lcl_input.GetNextSeqBatch(*lcl_scope));
#else
          CRef<CBlastQueryVector> query_batch(input.GetNextSeqBatch(*scope));
#endif


#if UPDATE
	    //#if USE_RESTART
	    //ceb increment completed batch counter for log file
	    batch_count++;//ceb
	    //#endif

	    // Get number of queries in batch
	    int lintNumberOfQueriesInBatch=query_batch->Size();

	    // Only team leader 0 needs to do this
	    if (local_tid==0) 
            {// If first batch

	      if (lintCurrentNumberOfThreadGroupQueries==0) {// Start at 0 
		lintpBatchStartPosition[thread_group_id] = 0;
	      } 
	      else{
		lintpBatchStartPosition[thread_group_id] = lintNextBatchStartPosition;
	      }
	      //if(tid==0)std::cout<<"checkpt 1.2 \n";std::cout.flush();
	      // Update running total 
	      lintCurrentNumberOfThreadGroupQueries += lintNumberOfQueriesInBatch;  
	      // Set query buffer end posiiton 
	      lintpBatchEndPosition[thread_group_id] = lintppThreadGroupQueryBufferEndPosition[thread_group_id][lintCurrentNumberOfThreadGroupQueries-1];
	      // Get ready for next iteraton 
	      lintNextBatchStartPosition = lintppThreadGroupQueryBufferEndPosition[thread_group_id][lintCurrentNumberOfThreadGroupQueries-1]+1;
	      //std::cout<<"checkpt 1.2.1 local_tid="<<local_tid<<" lintNextBatchStartPosition="<<lintNextBatchStartPosition<<"\n";std::cout.flush();
	    }
#endif

#if UPDATE
            #pragma omp barrier //wait for team leader thread to return window into query string (include conditional if(team leaders? > 0))
#endif //UPDATE
	  got_work = false; //reset flag
#if USE_RESTART
	  if ( !query_batch->Empty() && (batch_count > n_batches_precomputed) ) 
#else
	  if ( !query_batch->Empty() )
#endif//RESTART
	    {  got_work = true;}
	    
	  #pragma omp barrier //wait for shared variable got_work to be updated. Any thread returning true means true for all

	  CRef<CSearchResultSet> results;
          CLocalBlast* lcl_blast = NULL;

	  // No support for remote searches.
	  if ( got_work )  
          {
	    CRef<IQueryFactory> queries(new CObjMgr_QueryFactory(*query_batch));
	    SaveSearchStrategy(my_args, my_CmdLineArgs, queries, opts_hndl);

            //Create a local blast searcher for this thread
	    lcl_blast = new CLocalBlast(queries, opts_hndl, db_adapter);
	    lcl_blast->SetNumberOfThreads(my_CmdLineArgs->GetNumThreads());

#if UPDATE
	    lcl_blast->Run(); 
#else      
	    results = lcl_blast->Run();
#endif //UPDATE

#if UPDATE
               //Return objects saved from local_blast object
               //Collect pieces from local_blast object used to create output

	       //Make duplicates of these objects so we can delete the local_blast object
	       tHSPResults[tid] = gbhrpfncDuplicateBlastHSPResults(lcl_blast->GetHSPResults());
	       tQueryInfo[tid] = BlastQueryInfoDup(lcl_blast->GetQueryInfo());
	       tBlastScoreBlk[tid] = gbsbpfncDuplicateBlastScoreBlk(lcl_blast->GetScoreBlk());
  #if CLEAN_UP_OBJECTS
	       delete(lcl_blast);
  #endif //CLEAN

#endif //UPDATE

	  }
#if UPDATE
	  else
	    {
	      //Fill in non working thread contributions with some empty (but not NULL) structs
	      tHSPResults[tid] =    Blast_HSPResultsNew(0);
	      //results of BlastQueryInfoNew is not handled appropriately in combine function
	      //tQueryInfo[tid] =     BlastQueryInfoNew(eBlastTypeBlastx,0);//{0,0,0,NULL,0,NULL};
	      //tQueryInfo[tid] =     &BlastQueryInfo{0,0,0,NULL,0,NULL};
	      tQueryInfo[tid] =     new BlastQueryInfo{0,0,0,NULL,0,NULL};
	      tBlastScoreBlk[tid] = BlastScoreBlkNew(BLASTAA_SEQ_CODE,0);//type could be anything here
              query_batch->clear();
	    }

            #pragma omp barrier //ceb needed to merge results within thread group
	    // the results in parallel within the thread group.
	    if ( local_tid == 0 ){  
              //HSPResultsGroup[thread_group_id]=tHSPResults[tid];
	      //Concatenate and sort HSP results (by e-value) from threads in a thread group
	      HSPResultsGroup[thread_group_id] = gbhrpfncCombineTeamLeaderBlastHSPResults( tHSPResults,tid,num_team_leaders);
	      //These values should be the same inside a thread group so we only need to keep one thread
	      QueryInfoGroup[thread_group_id] = BlastQueryInfoDup(tQueryInfo[tid]);
	      //ceb ScoreBlk should be the same size but each corresponds to different database segment
	      //  so it's not certain how to merge these results for a thread group or what the ScoreBlk contains
	      BlastScoreBlkGroup[thread_group_id] = tBlastScoreBlk[tid];
	      //ceb a vector for consolidating query batches for output
	      //we only need one entry for each thread group
	      QueryBatchGroup[thread_group_id] = query_batch;
	    }
#endif

	  // Wait for all threads to finish. Then, reset all threads to 'see' the entire DB for reporting output.
          #pragma omp barrier//wait for all threads before changing db iteration range
	  my_seq_db->SetIterationRange(0,db_num_seqs);
          #pragma omp barrier //make sure all threads are here



#if UPDATE
  #if CLEAN_UP_OBJECTS
	    //free local thread objects
	    Blast_HSPResultsFree(tHSPResults[tid]);
	    if(tQueryInfo[tid]->num_queries>0){BlastQueryInfoFree(tQueryInfo[tid]);}//this function does not handle empty objects
  #endif
#endif //UPDATE


#if UPDATE
	    //we need to collect queries and store them in a big query_batch vector
            //This section should be called every time since at least one thread should have work here
	    if ( tid == 0 )
#endif //UPDATE
            { //Merge our consolidated thread group data across thread groups and output to file 
#if UPDATE
	      //thread info will be merged across thread groups prior to creating the results object           
	      BlastHSPResults* HSPResultsCollection = NULL; 
	      BlastQueryInfo* QueryInfoCollection = NULL;
	      BlastScoreBlk* ScoreBlkCollection = NULL; 
	      CRef<CBlastQueryVector> QueryBatchCollection(new CBlastQueryVector);  

	      HSPResultsCollection = gbhrpfncCombineThreadGroupBlastHSPResults(num_thread_groups,HSPResultsGroup);
	      QueryInfoCollection = gbqipfncCombineThreadGroupBlastQueryInfo(num_thread_groups,QueryInfoGroup);
	      ScoreBlkCollection = gbsbpfncCombineThreadGroupBlastScoreBlk(num_thread_groups,BlastScoreBlkGroup);
	      QueryBatchCollection = gbqvfncCombineReplicationGroupQueries(num_thread_groups,
									   buffer,
									   gbl_bsize,
									   lintpBatchStartPosition,
									   lintpBatchEndPosition,
									   &iconfig,
									   my_CmdLineArgs->GetQueryBatchSize(),
									   gbl_query_count,
									   gbl_scope);
#if USE_RESTART
//	       query_count += HSPResultsCollection->num_queries;//for restart
#endif
	       gbl_query_count += QueryBatchCollection->Size();
#endif //UPDATE

#if UPDATE
  #if CLEAN_UP_OBJECTS
	       //free objects merged within thread groups
	       for(int gid=0;gid<num_thread_groups;++gid)
	       {
		 Blast_HSPResultsFree(HSPResultsGroup[gid]);
		 if(QueryInfoGroup[gid]->num_queries>0){BlastQueryInfoFree(QueryInfoGroup[gid]);}
		 //BlastScoreBlkFree(BlastScoreBlkGroup[gid]);
		 //ok
		 //delete(BlastScoreBlkGroup[gid]);
	       }
  #endif //CLEAN
#endif //UPDATE

#if UPDATE
	       if(HSPResultsCollection->num_queries>0)
	       {
		 CRef<IQueryFactory> queries(new CObjMgr_QueryFactory(*QueryBatchCollection));
		 //This SeqInfoSrc needs to be created and filled with appropriate info
		 CRef<IBlastSeqInfoSrc> SeqInfoSrc;
		 SeqInfoSrc.Reset(db_adapter->MakeSeqInfoSrc());
       
		 results = gsrsfncCollectResults(HSPResultsCollection, 
						 QueryInfoCollection,
						 ScoreBlkCollection,
						 queries, 
						 &opt,
						 SeqInfoSrc);
		 BlastFormatter_PreFetchSequenceData(*results, gbl_scope);

		 //Loop over results
		 ITERATE(CSearchResultSet, result, *results)
		   {
  #if OUTPUT
		     formatter->PrintOneResultSet(**result, QueryBatchCollection);
#endif //OUTPUT
		   }
	       }
#else//UPDATE

	  // No support for archived writes.
	  if ( got_work )  {
	    BlastFormatter_PreFetchSequenceData(*results, scope);
	    ITERATE(CSearchResultSet, result, *results)
	      {
  #if OUTPUT
		  formatter->PrintOneResultSet(**result, query_batch);
  #endif//OUTPUT
	      }
	  }
#endif//UPDATE


	    
#if UPDATE
  #if CLEAN_UP_OBJECTS
	       //free objects merged across thread groups
	       //crashes code
	       //Blast_HSPResultsFree(HSPResultsCollection);
	       //delete(HSPResultsCollection);
	       //ok
	       BlastQueryInfoFree(QueryInfoCollection);
	       //crashes code
	       //BlastScoreBlkFree(ScoreBlkCollection);
  #endif //CLEAN
#endif //UPDATE
	    }//end if tid==0


          #pragma omp single
	       {working = false;}

	  #pragma omp critical
	  {
#if UPDATE
	    if ( !lcl_input.End() )
#else
	    if ( !input.End() )
#endif //UPDATE
	      {working = true;}
	  }
	  #pragma omp barrier

          #pragma omp flush(working)

	  #pragma omp barrier
        }//end working loop over superquerie
#if UPDATE
        if(tid==0)
#endif
//#if OUTPUT
	  {formatter->PrintEpilog(opt);}
//#endif
        if (my_CmdLineArgs->ProduceDebugOutput()) {
            opts_hndl->GetOptions().DebugDumpText(NcbiCerr, "BLAST options", 1);
        }

    } CATCH_ALL(status)

   }   // end the parallel block.
    return status;
}





#ifndef SKIP_DOXYGEN_PROCESSING
int main(int argc, char* argv[] /*, const char* envp[]*/)
{
  MPI_Init(&argc,&argv);

  MPI_Comm_rank(MPI_COMM_WORLD, &g_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &g_size);

  int host_name_len;
  char host_name[MPI_MAX_PROCESSOR_NAME];
  
  MPI_Get_processor_name( &(host_name[0]), &host_name_len );

  double start_time, end_time, runtime, maxruntime;

  start_time = MPI_Wtime();

  // Redirect output to a file.
  char new_out_stream_name[100];
  char new_err_stream_name[100];

  sprintf(new_out_stream_name, "hpc_blast.stdout.%05d", g_rank);
  sprintf(new_err_stream_name, "hpc_blast.stderr.%05d", g_rank);

  freopen( new_out_stream_name, "w", stdout );
  freopen( new_err_stream_name, "w", stderr );
  
  omp_set_nested(1);

  // Root rank (rank 0) does some checking on file existence and parameter correctness.
  int state = 1;    // Variable to indicate the state of the job. Assume its good until proven bad.
  int i = 0;        // Loop counter over argc.
  int j = 0;         // Additional loop counter.
  
  int db_index = 0; // The index in argv where the database name appears.
  int out_index = 0; //The index in argv where the output name is.

  char *db_name_for_rank = 0;
  char *output_file_name = 0;

  // Array for packing information to send all data in a single call.
  int bcast_data[5];

  if ( g_rank == 0 )
    {
      // Tasks:
      // 1) Parse out db name, query name, and num_threads parameters from command line.
      // 2) Check the existence and read the job parameters file.
      // 3) Sanity check on the values.
      // 4) Verify the existence of all files implied by the job parameters.
      // 5) Broadcast the job parameters to all ranks, or tell them to quit.

      // Task 1 : Parse out database name, query filename, and number of threads for the search engine from the command line.
      for ( i=0; i < argc; ++i )
	{
	  if ( (strcmp("-db",argv[i]) == 0 ) && (i+1) < argc )
	    {
	      db_str.assign(argv[i+1]);
	      db_index = i+1;
	      break;
	    }
	}

      if ( i == argc )  // didn't find a db option.
	{
	  std::cerr<<"FATAL ERROR: '-db' was not found on the command line."<<std::endl;
	  std::cerr.flush();
	  state = -1;
	}

      for ( i=0; i < argc; ++i )
	{
	  if ( (strcmp("-query",argv[i]) == 0 ) && (i+1) < argc )
	    {
	      query_str.assign(argv[i+1]);
	      break;
	    }
	}

      if ( i == argc )  // didn't find a query option.
	{
	  std::cerr<<"FATAL ERROR: '-query' was not found on the command line."<<std::endl;
	  std::cerr.flush();
	  state = -1;
	}

      for ( i=0; i < argc; ++i )
	{
	  if ( (strcmp("-num_threads",argv[i]) == 0 ) && (i+1) < argc )
	    {
	      num_threads_str.assign(argv[i+1]);
	      num_searcher_threads = atoi( argv[i+1] );
	      break;
	    }
	}

      if ( i == argc )  // didn't find a num_threads option.
	{
	  // Silence warnings.
	  //std::cerr<<"WARNING: '-num_threads' was not found on the command line. Assuming a value of 1."<<std::endl;
	  //std::cerr.flush();
	  
	  num_threads_str.assign("1");
	  num_searcher_threads = 1;
	}

      for ( i=0; i < argc; ++i )
	{
	  if ( (strcmp("-out",argv[i]) == 0 ) && (i+1) < argc )
	    {
	      out_str.assign(argv[i+1]);
	      out_index = i+1;
	      break;
	    }
	}

      if ( i == argc ) // didn't find the output file option
	{
	  std::cerr<<"FATAL ERROR: '-out' was not specified on the command line."<<std::endl;
	  std::cerr.flush();
	  state=-1;
	}

      // Task 2 : Check the existence of and read in the job parameters file.
      if ( file_exists("job_params") )
	{
	  std::ifstream job_input("job_params");
	  job_input >> num_thread_groups >> num_team_leaders >> num_ranks_per_group >> num_replication_groups;
	  job_input.close();
	}
      else
	{
	  state = -1;
	  std::cerr<<"FATAL ERROR: job_params file not found!"<<std::endl;
	  std::cerr.flush();
	}

      // Task 3 : Sanity check on the values.
      if ( state == 1 )
	{
	  // Check num_searcher_threads.
	  if ( num_searcher_threads < 1 )
	    {
	      std::cerr<<"FATAL ERROR: num_searcher_threads must be positive."<<std::endl;
	      std::cerr<<"  num_searcher_threads = "<<num_searcher_threads<<std::endl;
	      std::cerr.flush();
	      state = -1;
	    }
	  else if ( num_searcher_threads > 100 ) // Max db partitions by the search engine.
	    {
	      std::cerr<<"FATAL ERROR: The number of searcher threads for the search engine exceeds the number of"<<std::endl;
	      std::cerr<<"             database partitions the engine supplies."<<std::endl;
	      std::cerr.flush();
	      state = -1;
	    }

	  // Check num_thread_groups.
	  if ( num_thread_groups < 1 )
	    {
	      std::cerr<<"FATAL ERROR: num_thread_groups must be positive."<<std::endl;
	      std::cerr<<"  num_thread_groups = "<<num_thread_groups<<std::endl;
	      std::cerr.flush();
	      state = -1;
	    }
	  else if ( num_thread_groups > 240 && num_thread_groups < 361 )
	    {
	      std::cerr<<"WARNING: num_thread_groups will oversubscribe resources, but not too high."<<std::endl;
	      std::cerr<<"  240 < (num_thread_groups="<<num_thread_groups<<") <= 360"<<std::endl;
	      std::cerr.flush();
	    }
	  else if ( num_thread_groups > 360 )
	    {
	      std::cerr<<"FATAL ERROR: num_thread_groups is too large."<<std::endl;
	      std::cerr<<"  360 < num_thread_groups="<<num_thread_groups<<std::endl;
	      std::cerr.flush();
	      state = -1;
	    }

	  // Check num_team_leaders.
	  if ( num_team_leaders < 1 )
	    {
	      std::cerr<<"FATAL ERROR: num_team_leaders must be positive."<<std::endl;
	      std::cerr<<"  num_team_leaders = "<<num_team_leaders<<std::endl;
	      std::cerr.flush();
	      state = -1;
	    }
	  else if ( num_team_leaders > 240 && num_team_leaders < 361 )
	    {
	      std::cerr<<"WARNING: num_team_leaders will oversubscribe resources, but not too high."<<std::endl;
	      std::cerr<<"  240 < (num_team_leaders="<<num_team_leaders<<") <= 360"<<std::endl;
	      std::cerr.flush();
	    }
	  else if ( num_team_leaders > 360 )
	    {
	      std::cerr<<"FATAL ERROR: num_team_leaders is too large."<<std::endl;
	      std::cerr<<"  360 < num_team_leaders="<<num_team_leaders<<std::endl;
	      std::cerr.flush();
	      state = -1;
	    }
	  
	  // Check the total number of threads used by a rank.
	  num_threads_total = num_thread_groups * num_team_leaders * num_searcher_threads ;
	  if ( num_threads_total < 1 )
	    {
	      std::cerr<<"FATAL ERROR: num_threads_total must be positive."<<std::endl;
	      std::cerr<<"  num_threads_total = "<<num_threads_total<<std::endl;
	      std::cerr.flush();
	      state = -1;
	    }
	  else if ( num_threads_total > 240 && num_threads_total < 361 )
	    {
	      std::cerr<<"WARNING: num_threads_total will oversubscribe resources, but not too high."<<std::endl;
	      std::cerr<<"  240 < (num_threads_total="<<num_threads_total<<") <= 360"<<std::endl;
	      std::cerr.flush();
	    }
	  else if ( num_threads_total > 360 )
	    {
	      std::cerr<<"FATAL ERROR: num_threads_total is too large."<<std::endl;
	      std::cerr<<"  360 < num_threads_total="<<num_threads_total<<std::endl;
	      std::cerr.flush();
	      state = -1;
	    }

	  // Test MPI values.
	  if ( num_replication_groups < 1 )
	    {
	      std::cerr<<"FATAL ERROR: num_replication_groups must be positive."<<std::endl;
	      std::cerr<<"  num_replication_groups = "<<num_replication_groups<<std::endl;
	      std::cerr.flush();
	      state = -1;
	    }

	  if ( num_ranks_per_group < 1 )
	    {
	      std::cerr<<"FATAL ERROR: num_ranks_per_group must be positive."<<std::endl;
	      std::cerr<<"  num_ranks_per_group = "<<num_ranks_per_group<<std::endl;
	      std::cerr.flush();
	      state = -1;
	    }

	  num_ranks_total = num_replication_groups * num_ranks_per_group;

	  if ( num_ranks_total != g_size )
	    {
	      std::cerr<<"FATAL ERROR: num_ranks_total requested by the job does not match number of ranks in COMM_WORLD."<<std::endl;
	      std::cerr<<"  num_ranks_total = "<<num_ranks_total<<std::endl;
	      std::cerr<<"  Size of MPI_COMM_WORLD = "<<g_size<<std::endl;
	      std::cerr.flush();
	      state = -1;
	    }

	} // End if valid state value for Task 3.

      // Task 4 : Verify existence of all files that the parameters imply should exist.
      //if ( state == 1 )
      if ( 0 )
	{
	  // First check database files.
	  std::string file_name;

	  if ( num_ranks_per_group == 1 )
	    {
	      file_name.assign(db_str);

	      file_name += ".pin";
	      if ( ! file_exists(file_name) )
		{
		  std::cerr<<"FATAL ERROR: The database file does not exist: "<<file_name<<std::endl;
		  std::cerr.flush();
		  state = -1;
		}

	      file_name.assign(db_str);
	      file_name += ".psq";
	      if ( ! file_exists(file_name) )
		{
		  std::cerr<<"FATAL ERROR: The database file does not exist: "<<file_name<<std::endl;
		  std::cerr.flush();
		  state = -1;
		}

	      file_name.assign(db_str);
	      file_name += ".phr";
	      if ( ! file_exists(file_name) )
		{
		  std::cerr<<"FATAL ERROR: The database file does not exist: "<<file_name<<std::endl;
		  std::cerr.flush();
		  state = -1;
		}
	    }
	  else  // multiple partitions.
	    {
	      for ( i=0; i < num_ranks_per_group; ++i )
		{
		  std::stringstream sstm;

		  sstm << db_str;

		  if ( i < 10 )
		    sstm << ".0" << i;
		  else
		    sstm << "." << i;

		  file_name.assign( sstm.str() + ".pin" );
		  if ( ! file_exists(file_name) )
		    {
		      std::cerr<<"FATAL ERROR: The database file does not exist: "<<file_name<<std::endl;
		      std::cerr.flush();
		      state = -1;
		    }

		  file_name.assign( sstm.str() + ".psq" );
		  if ( ! file_exists(file_name) )
		    {
		      std::cerr<<"FATAL ERROR: The database file does not exist: "<<file_name<<std::endl;
		      std::cerr.flush();
		      state = -1;
		    }

		  file_name.assign( sstm.str() + ".phr" );
		  if ( ! file_exists(file_name) )
		    {
		      std::cerr<<"FATAL ERROR: The database file does not exist: "<<file_name<<std::endl;
		      std::cerr.flush();
		      state = -1;
		    }
		}
	    }

	  // Now check query file. Assume we will always partition the queries ( or at least pack them so they will always be in the given file naming convention.

	  for ( i=0; i < num_replication_groups; ++i )
	    {
	      
	      for ( j=0; j < num_thread_groups; ++j )
		{
		  std::stringstream sstm;

		  sstm << query_str << "." << i << "." << j;

		  file_name = sstm.str();

		  if ( ! file_exists(file_name) )
		    {
		      std::cerr<<"FATAL ERROR: The query file does not exist: "<<file_name<<std::endl;
		      std::cerr.flush();
		      state = -1;
		    }
		}
	      
	    }
	  
	} // End if valid state value for Task 4.

    } // End if root rank.

  // Task 5 : Broadcast the state of the job to all ranks. If a valid job, broadcast parameters to all ranks.

  if ( g_rank == 0 ) 
    {
      std::cout<<"Root is about to broadcast the state variable: "<<state<<std::endl;
      std::cout.flush();
    }
  
  MPI_Bcast ( &state, 1, MPI_INT, 0, MPI_COMM_WORLD );

  if ( state == -1 )
    {
      MPI_Finalize();
      return 1;
    }
  else
    {
      // Root packs the data for broadcast.
      if ( g_rank == 0 )
	{
	  bcast_data[0] = num_thread_groups;
	  bcast_data[1] = num_team_leaders;
	  bcast_data[2] = num_searcher_threads;
	  bcast_data[3] = num_ranks_per_group;
	  bcast_data[4] = num_replication_groups;

	  std::cout<<"Sending bcast_data[0] = num_thread_groups      = "<<num_thread_groups<<std::endl;
	  std::cout<<"Sending bcast_data[0] = num_team_leaders       = "<<num_team_leaders<<std::endl;
	  std::cout<<"Sending bcast_data[0] = num_searcher_threads   = "<<num_searcher_threads<<std::endl;
	  std::cout<<"Sending bcast_data[0] = num_ranks_per_group    = "<<num_ranks_per_group<<std::endl;
	  std::cout<<"Sending bcast_data[0] = num_replication_groups = "<<num_replication_groups<<std::endl;
	  std::cout.flush();
	}

      MPI_Bcast ( bcast_data, 5, MPI_INT, 0, MPI_COMM_WORLD );

      if ( g_rank != 0 )
	{
	  // Unpack the sent data.
	  num_thread_groups      = bcast_data[0];
	  num_team_leaders       = bcast_data[1];
	  num_searcher_threads   = bcast_data[2];
	  num_ranks_per_group    = bcast_data[3];
	  num_replication_groups = bcast_data[4];

	  std::cout<<"Received bcast_data[0] = num_thread_groups      = "<<num_thread_groups<<std::endl;
	  std::cout<<"Received bcast_data[0] = num_team_leaders       = "<<num_team_leaders<<std::endl;
	  std::cout<<"Received bcast_data[0] = num_searcher_threads   = "<<num_searcher_threads<<std::endl;
	  std::cout<<"Received bcast_data[0] = num_ranks_per_group    = "<<num_ranks_per_group<<std::endl;
	  std::cout<<"Received bcast_data[0] = num_replication_groups = "<<num_replication_groups<<std::endl;
	  std::cout.flush();

	  // Also get the database name and query file name from the command line. We know it is there as implied by the state.
	  for ( i=0; i < argc; ++i )
	    {
	      if ( (strcmp("-db",argv[i]) == 0 ) && (i+1) < argc )
		{
		  db_str.assign(argv[i+1]);
		  db_index = i+1;
		  break;
		}
	    }

	  for ( i=0; i < argc; ++i )
	    {
	      if ( (strcmp("-query",argv[i]) == 0 ) && (i+1) < argc )
		{
		  query_str.assign(argv[i+1]);
		  break;
		}
	    }

	  for ( i=0; i < argc; ++i )
	    {
	      if ( (strcmp("-out",argv[i]) == 0 ) && (i+1) < argc )
		{
		  out_str.assign(argv[i+1]);
		  out_index= i+1;
		  break;
		}
	    }
	  
	}
    }

  MPI_Barrier( MPI_COMM_WORLD );

  // Establish the rank's ID in the replication group.
  rep_rank = g_rank % num_ranks_per_group;          // Which db partition to read.
  rep_group = g_rank / num_ranks_per_group;         // Which query partition to read.

  // If the database is partitioned, we will have to change the process' db file name so that it will open the appropriate file.
  if ( num_ranks_per_group > 1 )
    {
      db_name_for_rank = new char[ strlen(argv[db_index]) + 10 ];

      if ( rep_rank < 10 )
	sprintf(db_name_for_rank,"%s.0%d",argv[db_index],rep_rank);
      else
	sprintf(db_name_for_rank,"%s.%d",argv[db_index],rep_rank);

      argv[db_index] = db_name_for_rank;
    }

  // Allow for each rank to have its own output file. This is purely for the sake of the NCBI toolkit which will attempt to open it for writing later (which is overridden in the threaded section).
  if ( g_size > 1 )
    {
      output_file_name = new char[ strlen(argv[out_index]) + 10];

      sprintf(output_file_name,"%s.%d",argv[out_index],g_rank);

      argv[out_index] = output_file_name;
    }

  std::cout<<"Rank "<<g_rank<<" is running on "<<host_name<<std::endl;
  std::cout<<"Rank "<<g_rank<<" : my database is "<<argv[db_index]<<std::endl;

  // Start the BLAST search.

  int blast_return = CBlastxApp().AppMain(argc, argv, 0, eDS_Default, 0);

  end_time = MPI_Wtime();

  runtime = end_time - start_time;

  std::cout<<"Runtime is "<<runtime<<std::endl;
  std::cout.flush();

  // Catch all nodes at the end.
  MPI_Barrier(MPI_COMM_WORLD);

  MPI_Reduce( &runtime, &maxruntime, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD );

  if ( g_rank == 0 )
    std::cout<<"Job runtime is "<<maxruntime<<std::endl;

  MPI_Finalize();

  if ( db_name_for_rank != 0 )
    delete [] db_name_for_rank;

  if ( output_file_name != 0 )
    delete [] output_file_name;

  return blast_return;
}











// Functions to combine objects across within thread groups ==================================

BlastHSPResults* gbhrpfncCombineTeamLeaderBlastHSPResults ( char ** HSPResultsBuffers,
							    size_t * buffer_size,
							    int tid,
							    int num_team_leaders )
{
  Int4 num_queries = 0;
  Int4 dummy_int4 = 0;  // dummy variable.
  double dummy_double = 0.;
  Boolean dummy_boolean;
  Boolean all_heapified = 0;

  Int4 total_hsplist_current = 0;
  Int4 total_hsplist_count = 0;
  Int4 total_hsplist_max=0;
  Int4 total_hsparray_allocated = 0;

  double worst_evalue_all=0.;
  Int4 low_score_all=0;

  // Need a buffer pointer for each different buffer.
  char **buff_ptr = NULL;

  Int4 *hsplist_counts = NULL;
  //ceb
  Int4 *total_hsplist_counts = NULL;

  buff_ptr = (char**)malloc(num_team_leaders * sizeof(char*));
  hsplist_counts = (Int4*)malloc(num_team_leaders * sizeof(char*));

  //std::cout<<std::endl;//ceb

  BlastHSPResults *hsp_results = NULL;    // The guy we're returning.

  // Before we unpack the buffers into the object we have to know how many of certain objects their are.
  // So, first we crawl across the buffers to count the stuff up.
  for ( int i=0; i < num_team_leaders; ++i )
    {
      //std::cout<< "chkpt -1 tid: " << tid <<" i: "<< i <<std::endl;std::cout.flush();
      buff_ptr[i] = HSPResultsBuffers[tid+i];
    }


  // Get the number of queries -- its the same across all threads in a thread group.
  memcpy( &num_queries , buff_ptr[0], sizeof(Int4) );
  //ceb Seems to fix the offset mismatch error when going in to query 0+1
  for ( int i=0; i < num_team_leaders; ++i )
    { buff_ptr[i] += sizeof(Int4); }//each thread will have num_queries, we only read it in once outside this loop and then skip here.


  //ceb
  //total_hsplist_counts = (Int4*)malloc(num_queries * sizeof(char*));
  //for(int i=0;i<num_queries;++i){total_hsplist_counts[i]=0;}


  //if(tid>=0){std::cout<<"tid: "<<tid<<" Found "<<num_queries<<" num_queries in the first buffer."<<std::endl;}

  hsp_results = Blast_HSPResultsNew(num_queries);

  // We have to move across the buffer in a very particular way (reverse of how its packed.)
  // See traceback_stage.cpp for details on how its packed.

  for ( int BlastHSPResults_query_index = 0; BlastHSPResults_query_index < hsp_results->num_queries ; ++BlastHSPResults_query_index )
    {    
      // Allocate the BlastHitList for this query.
      BlastHitList *bhl = Blast_HitListNew(1);

      //if(tid>=0){std::cout<<"tid: "<<tid<<" Query "<<BlastHSPResults_query_index<<std::endl;std::cout.flush();}

      // hsplist_count
      total_hsplist_count=0;//ceb Reset here for each query
      for ( int i=0; i < num_team_leaders; ++i )
	{
	  // read past num_queries
	  //buff_ptr[i] += sizeof(Int4);//each thread will have num_queries, we only read it in once outside this loop and then skip here.
          memcpy( &dummy_int4, buff_ptr[i], sizeof(Int4) );
	  buff_ptr[i] += sizeof(Int4);
	  
	  //if(tid>=0){std::cout<<"  "<<"Team leader "<<i<<" has hsplist_count "<<dummy_int4<<std::endl; std::cout.flush();}

          //ceb We want to merge hsplist data across multiple team leaders in a thread group for a given query.
	  total_hsplist_count += dummy_int4;

	  //total_hsplist_count[BlastHSPResults_query_index] += dummy_int4;//ceb

	  hsplist_counts[i] = dummy_int4;
	}

      //if(tid>=0){std::cout<<"tid: "<<tid<<" Query "<<BlastHSPResults_query_index<<" Total hsplist_count "<<total_hsplist_count<<std::endl; std::cout.flush();}

      bhl->hsplist_count = total_hsplist_count;

      //if(tid==0){std::cout<<"checkpt 1 "<<std::endl; std::cout.flush();}

      // hsplist_max
      for ( int i=0; i < num_team_leaders; ++i )
	{
	  memcpy( &dummy_int4, buff_ptr[i], sizeof(Int4) );
	  buff_ptr[i] += sizeof(Int4);

	  //if(tid>=0){std::cout<<"  "<<"Team leader "<<i<<" has hsplist_max "<<dummy_int4<<std::endl; std::cout.flush();}

	  if ( i == 0 )
	    {
	      //memcpy( &dummy_int4, buff_ptr[i], sizeof(Int4) );
	      total_hsplist_max = dummy_int4;
	    }
	}

      //if(tid>=0){std::cout<<"tid: "<<tid<<" pre mult total_hsplist_max: "<<total_hsplist_max<< std::endl; std::cout.flush();}

      while ( total_hsplist_max < total_hsplist_count )
	{total_hsplist_max *= 2;}

      //if(tid>=0){std::cout<<"tid: "<<tid<<" post mult Total hsplist max is "<<total_hsplist_max<<std::endl;std::cout.flush();}

      bhl->hsplist_max = total_hsplist_max;

      // worst evalue
      memcpy( &dummy_double, buff_ptr[0], sizeof(double));
      buff_ptr[0] += sizeof(double);
      worst_evalue_all = dummy_double;
      //if(tid>=0){std::cout<<"tid: "<<tid<<"checkpt 4 "<<std::endl; std::cout.flush();}
      //if(tid>=0){std::cout<<"tid: "<<tid<<"  "<<"Team leader "<<"0"<<" has worst evalue "<<dummy_double<<std::endl;std::cout.flush();}
      for ( int i=1; i < num_team_leaders; ++i )
	{
	  memcpy( &dummy_double, buff_ptr[i], sizeof(double) );
	  buff_ptr[i] += sizeof(double);
	  //if(tid>=0){std::cout<<"tid: "<<tid<<"  "<<"Team leader "<<i<<" has worst evalue "<<dummy_double<<std::endl;std::cout.flush();}
	  if( dummy_double > worst_evalue_all )
	    {worst_evalue_all = dummy_double;}
	}
      //if(tid>=0){std::cout<<"tid: "<<tid<<"  "<<"Worst evalue all is "<<worst_evalue_all<<std::endl;std::cout.flush();}
      bhl->worst_evalue = worst_evalue_all;


      // best bit(low) score
      memcpy( &dummy_int4, buff_ptr[0], sizeof(Int4));
      buff_ptr[0] += sizeof(Int4);
      low_score_all = dummy_int4;
      //if(tid>=0){std::cout<<"tid: "<<tid<<"checkpt 6 "<<std::endl; std::cout.flush();}
      //if(tid>=0){std::cout<<"  "<<"Team leader "<<"0"<<" has low score "<<dummy_int4<<std::endl;std::cout.flush();}
      for ( int i=1; i < num_team_leaders; ++i )
	{
	  memcpy( &dummy_int4, buff_ptr[i], sizeof(Int4));
	  buff_ptr[i] += sizeof(Int4);
	  //if(tid>=0){std::cout<<"tid: "<<tid<<"  "<<"Team leader "<<i<<" has low score "<<dummy_int4<<std::endl;std::cout.flush();}
	  if( dummy_int4 < low_score_all )
	    {low_score_all = dummy_int4;}
	}
      //if(tid>=0){std::cout<<"tid: "<<tid<<"  "<<"Low score all is "<<low_score_all<<std::endl;std::cout.flush();}
      bhl->low_score = low_score_all;

      //heapified boolean
      for ( int i=0; i < num_team_leaders; ++i )
	{
	  memcpy( &dummy_boolean, buff_ptr[i], sizeof(Boolean) );
	  buff_ptr[i] += sizeof(Boolean);
	  //if(tid>=0)std::cout<<"tid: "<<tid<<"  "<<"Team leader "<<i<<" has heapified "<<(int)dummy_boolean<<std::endl;std::cout.flush();
	  if ( dummy_boolean == 1 )
	    {all_heapified = dummy_boolean;}
	}
      // if any are heapified assume the merge can be.
      //if(tid>=0){std::cout<<"tid: "<<tid<<"  "<<"all heapified is "<<(int)all_heapified<<std::endl;std::cout.flush();}
      bhl->heapified = all_heapified;

      //hsplist_current
      for ( int i=0; i < num_team_leaders; ++i )
	{
	  memcpy( &dummy_int4, buff_ptr[i], sizeof(Int4) );
	  buff_ptr[i] += sizeof(Int4);
	  //if(tid>=0){std::cout<<"tid: "<<tid<<"  "<<"Team leader "<<i<<" has hsplist_current "<<dummy_int4<<std::endl;std::cout.flush();}
	}
      //if(tid>=0){std::cout<<"tid: "<<tid<<"checkpt 10 "<<std::endl; std::cout.flush();}


      // take the first one or the max. add 100 until it surpasses hsplist_count!
      while ( dummy_int4 < total_hsplist_count )
	{dummy_int4 += 100;}

      total_hsplist_current = dummy_int4;
      //if(tid>=0){std::cout<<"tid: "<<tid<<"  Total hsplist_current is "<<total_hsplist_current<<std::endl; std::cout.flush();}
      bhl->hsplist_current = total_hsplist_current;

      // Allocate the hspList to current long.
      bhl->hsplist_array = (BlastHSPList**)calloc( bhl->hsplist_current, sizeof(BlastHSPList*) );

      //if(tid>=0){std::cout<<"tid: "<<tid<<"checkpt 11 "<<std::endl; std::cout.flush();}

      Int4 global_hsplist_index = 0;
      // Loop over each team leader and process the hit lists.
      for ( int i=0; i < num_team_leaders; ++i )
	{
	  //if(tid>=0){std::cout<<"tid: "<<tid<<"checkpt 12 hsplist_counts["<<i<<"]: "<<hsplist_counts[i]<<std::endl; std::cout.flush();}

	  //for ( Int4 BlastHitList_hsplist_index; BlastHitList_hsplist_index < hsplist_counts[i]; ++BlastHitList_hsplist_index )
	  for ( int BlastHitList_hsplist_index=0; BlastHitList_hsplist_index < hsplist_counts[i]; ++BlastHitList_hsplist_index )
	    {
	      // Allocate a BlastHSPList ptr
	      BlastHSPList* bhspl = (BlastHSPList*) calloc(1,sizeof(BlastHSPList) );

	      memcpy( &(bhspl->oid) , buff_ptr[i], sizeof(Int4) );
	      buff_ptr[i] += sizeof(Int4);
	      //if(tid>=0){std::cout<<"tid: "<<tid<<" bhspl->oid: "<<bhspl->oid<<std::endl; std::cout.flush();}
	      
	      memcpy( &(bhspl->query_index) , buff_ptr[i], sizeof(Int4) );
	      buff_ptr[i] += sizeof(Int4);

	      memcpy( &(bhspl->hspcnt) , buff_ptr[i], sizeof(Int4) );
	      buff_ptr[i] += sizeof(Int4);

	      memcpy( &(bhspl->allocated) , buff_ptr[i], sizeof(Int4) );
	      buff_ptr[i] += sizeof(Int4);

	      memcpy( &(bhspl->hsp_max) , buff_ptr[i], sizeof(Int4) );
	      buff_ptr[i] += sizeof(Int4);
	      //if(tid>=0){std::cout<<"tid: "<<tid<<" bhspl->hsp_max: "<<bhspl->hsp_max<<std::endl; std::cout.flush();}

	      memcpy( &(bhspl->do_not_reallocate) , buff_ptr[i], sizeof(Boolean) );
	      buff_ptr[i] += sizeof(Boolean);

	      memcpy( &(bhspl->best_evalue) , buff_ptr[i], sizeof(double) );
	      buff_ptr[i] += sizeof(double);
	      //if(tid>=0){std::cout<<"tid: "<<tid<<" bhspl->best_evalue "<<bhspl->best_evalue<<std::endl; std::cout.flush();}

	      // Allocate the array for High Scoring segment pairs.
	      bhspl->hsp_array = (BlastHSP**) calloc( bhspl->allocated, sizeof(BlastHSP*) );
	      //if(tid>=0){std::cout<<"tid: "<<tid<<"bhspl->hspcnt "<<bhspl->hspcnt<<std::endl; std::cout.flush();}

	      // Loop over all HSPs between the query and the subject sequence.
	      for ( int BlastHSP_hsp_array_index=0; BlastHSP_hsp_array_index < bhspl->hspcnt; ++BlastHSP_hsp_array_index )
		{
		  BlastHSP* bhsp = (BlastHSP*) calloc (1, sizeof(BlastHSP) );

		  memcpy ( &(bhsp->score), buff_ptr[i] , sizeof(Int4) );
		  buff_ptr[i] += sizeof(Int4);
		  
		  memcpy ( &(bhsp->num_ident), buff_ptr[i] , sizeof(Int4) );
		  buff_ptr[i] += sizeof(Int4);
		
		  memcpy ( &(bhsp->bit_score), buff_ptr[i] , sizeof(double) );
		  buff_ptr[i] += sizeof(double);

		  memcpy ( &(bhsp->evalue), buff_ptr[i] , sizeof(double) );
		  buff_ptr[i] += sizeof(double);

		  memcpy ( &(bhsp->context), buff_ptr[i] , sizeof(Int4) );
		  buff_ptr[i] += sizeof(Int4);

		  memcpy ( &(bhsp->num), buff_ptr[i] , sizeof(Int4) );
		  buff_ptr[i] += sizeof(Int4);

		  memcpy ( &(bhsp->comp_adjustment_method), buff_ptr[i] , sizeof(Int2) );
		  buff_ptr[i] += sizeof(Int2);

		  memcpy ( &(bhsp->num_positives), buff_ptr[i] , sizeof(Int4) );
		  buff_ptr[i] += sizeof(Int4);

		  memcpy ( &(bhsp->query.frame), buff_ptr[i] , sizeof(Int2) );
		  buff_ptr[i] += sizeof(Int2);

		  memcpy ( &(bhsp->query.offset), buff_ptr[i] , sizeof(Int4) );
		  buff_ptr[i] += sizeof(Int4);

		  memcpy ( &(bhsp->query.end), buff_ptr[i] , sizeof(Int4) );
		  buff_ptr[i] += sizeof(Int4);

		  memcpy ( &(bhsp->query.gapped_start), buff_ptr[i] , sizeof(Int4) );
		  buff_ptr[i] += sizeof(Int4);

		  memcpy ( &(bhsp->subject.frame), buff_ptr[i] , sizeof(Int2) );
		  buff_ptr[i] += sizeof(Int2);

		  memcpy ( &(bhsp->subject.offset), buff_ptr[i] , sizeof(Int4) );
		  buff_ptr[i] += sizeof(Int4);

		  memcpy ( &(bhsp->subject.end), buff_ptr[i] , sizeof(Int4) );
		  buff_ptr[i] += sizeof(Int4);

		  memcpy ( &(bhsp->subject.gapped_start), buff_ptr[i] , sizeof(Int4) );
		  buff_ptr[i] += sizeof(Int4);
		  //if(tid>=0){std::cout<<"bhsp->subject.gapped_start "<<bhsp->subject.gapped_start<<std::endl; std::cout.flush();}

		  // The next value might be a setinel value.
		  Int4 sentinel;
                  
 
		  memcpy ( &sentinel , buff_ptr[i], sizeof(Int4) );

		  if ( sentinel )
		    {
		      SPHIHspInfo *my_pat_info = (SPHIHspInfo*) calloc( 1, sizeof(SPHIHspInfo) );

		      memcpy( &(my_pat_info->index) , buff_ptr[i], sizeof(Int4) );
		      buff_ptr[i] += sizeof(Int4);

		      memcpy( &(my_pat_info->length) , buff_ptr[i], sizeof(Int4) );
		      //if(tid>=0){std::cout<<"tid: "<<tid<<"my_pat_info->length "<<my_pat_info->length<<std::endl; std::cout.flush();}
		      buff_ptr[i] += sizeof(Int4);

		      bhsp->pat_info = my_pat_info ;
		    }
		  else
		    {
		      bhsp->pat_info = NULL;
		      buff_ptr[i] += sizeof(Int4);
		    }

		  // The next value might be a sentinel value on GapEditScript.
		  memcpy ( &sentinel , buff_ptr[i] , sizeof(Int4) );

		  if ( sentinel )
		    {
		      GapEditScript *my_gap_info = (GapEditScript*) calloc(1, sizeof(GapEditScript) );

		      memcpy( &(my_gap_info->size) , buff_ptr[i] , sizeof(Int4) );
		      //if(tid>=0){std::cout<<"tid: "<<tid<<" my_gap_info->size "<<my_gap_info->size<<std::endl; std::cout.flush();}

		      buff_ptr[i] += sizeof(Int4);

		      my_gap_info->op_type = (EGapAlignOpType*) calloc( my_gap_info->size , sizeof(EGapAlignOpType) );
		      memcpy( &(my_gap_info->op_type[0]), buff_ptr[i] , my_gap_info->size * sizeof(EGapAlignOpType) );
		      buff_ptr[i] += ( my_gap_info->size * sizeof(EGapAlignOpType) );

		      my_gap_info->num = (Int4*) calloc( my_gap_info->size , sizeof(Int4) );
		      memcpy( &(my_gap_info->num[0]), buff_ptr[i] , my_gap_info->size * sizeof(Int4) );
		      buff_ptr[i] += (my_gap_info->size * sizeof(Int4) );

		      bhsp->gap_info = my_gap_info ;
		    }
		  else
		    {
		      bhsp->gap_info = NULL; 
		      buff_ptr[i] += sizeof(Int4);
		    }

		  // Now attach this BlastHSP to the parent BlastHSPList.
		  bhspl->hsp_array[ BlastHSP_hsp_array_index ] = bhsp ;

		} // End of loop over the BlastHSP_hsp_array_index -> the high scoring pair data

              //if(tid==0){std::cout<<"tid: "<<tid<<" checkpt 19 "<<std::endl; std::cout.flush();}

	      // Now attach this BlastHSPList to the parent BlastHitList
	      bhl->hsplist_array[ global_hsplist_index ] = bhspl ;

	      ++global_hsplist_index;  // increment the global index
	      
	    }  // End loop over BlastHitList_hsplist_index

	} // End Loop over team leaders.

      //if(tid==0){std::cout<<"tid: "<<tid<<" checkpt 20 "<<std::endl; std::cout.flush();}

      // Now we can add this BlastHitList to the parent BlastHSPResults
      hsp_results->hitlist_array[ BlastHSPResults_query_index ] = bhl ;

  }

  //if(tid>=0){std::cout<<"tid: "<<tid<<" final hsp_results->num_queries: "<<hsp_results->num_queries <<std::endl; std::cout.flush();}
  // for ( int BlastHSPResults_query_index = 0; BlastHSPResults_query_index < hsp_results->num_queries ; ++BlastHSPResults_query_index )
  //  {
  //    if(tid>=0){std::cout<<"tid: "<<tid<<" Query "<<BlastHSPResults_query_index<<std::endl; std::cout.flush();}
  //    //if(tid>=0){std::cout<<"tid: "<<tid<<" final hsp_results->hitlist_array["<< BlastHSPResults_query_index <<"]->worst_evalue: "<< hsp_results->hitlist_array[BlastHSPResults_query_index]->worst_evalue <<std::endl; std::cout.flush();}
  //    if(tid>=0){std::cout<<"tid: "<<tid<<" final hsp_results->hitlist_array["<< BlastHSPResults_query_index <<"]->hsplist_count: "<< hsp_results->hitlist_array[BlastHSPResults_query_index]->hsplist_count <<std::endl; std::cout.flush();}
  //    //if(tid>=0){std::cout<<"tid: "<<tid<<" final hsp_results->hitlist_array["<< BlastHSPResults_query_index <<"]->hsplist_max: "<< hsp_results->hitlist_array[BlastHSPResults_query_index]->hsplist_max <<std::endl; std::cout.flush();}
  //    //if(tid>=0){std::cout<<"tid: "<<tid<<" final hsp_results->hitlist_array["<< BlastHSPResults_query_index <<"]->low_score: "<< hsp_results->hitlist_array[BlastHSPResults_query_index]->low_score <<std::endl; std::cout.flush();}
  //    //if(tid>=0){std::cout<<"tid: "<<tid<<" final hsp_results->hitlist_array["<< BlastHSPResults_query_index <<"]->heapified: "<< hsp_results->hitlist_array[BlastHSPResults_query_index]->heapified <<std::endl; std::cout.flush();}
  //    if(tid>=0){std::cout<<"tid: "<<tid<<" final hsp_results->hitlist_array["<< BlastHSPResults_query_index <<"]->hsplist_current: "<< hsp_results->hitlist_array[BlastHSPResults_query_index]->hsplist_current <<std::endl; std::cout.flush();}
  //  }



  //if(tid>=0){std::cout<<"tid: "<<tid<<" Pre SortByValue "<<std::endl; std::cout.flush();}
  // Now we can sort the HSPResults.
  Blast_HSPResultsSortByEvalue( hsp_results ) ;
  //if(tid>=0){std::cout<<"tid: "<<tid<<" Post SortByValue "<<std::endl; std::cout.flush();}


  return hsp_results;

}








BlastHSPResults* gbhrpfncCombineTeamLeaderBlastHSPResults ( BlastHSPResults** HSPResultsArray,
							    int tid,
							    int num_team_leaders )
{
  Int4 num_queries = 0;
  Int4 dummy_int4 = 0;  // dummy variable.
  double dummy_double = 0.;
  Boolean dummy_boolean;
  Boolean all_heapified = 0;

  Int4 total_hsplist_current = 0;
  Int4 total_hsplist_count = 0;
  Int4 total_hsplist_max=0;
  Int4 total_hsparray_allocated = 0;

  double worst_evalue_all=0.;
  Int4 low_score_all=0;

  // Need a buffer pointer for each different buffer.
  Int4 *hsplist_counts = NULL;
  BlastHSPResults **hsp_results =  new BlastHSPResults*[num_team_leaders];    // use this to point to a value in HSPResultsArray
  BlastHSPResults *merged_hsp_results = NULL;    // The guy we're returning.

  hsplist_counts = (Int4*)malloc(num_team_leaders * sizeof(char*));

  for ( int i=0; i < num_team_leaders; ++i )
  {
    hsp_results[i] = HSPResultsArray[tid+i];
  }

  // Get the number of queries -- its the same across all threads in a thread group.
  num_queries = (hsp_results[0])->num_queries;
  //std::cout<<"(hsp_results[0])->num_queries= "<<num_queries <<std::endl;std::cout.flush();
  merged_hsp_results = Blast_HSPResultsNew(num_queries);

  // We have to move across the buffer in a very particular way (reverse of how its packed.)
  // See traceback_stage.cpp for details on how its packed.

  for ( int BlastHSPResults_query_index = 0; BlastHSPResults_query_index < merged_hsp_results->num_queries ; ++BlastHSPResults_query_index )
    {    
      // Allocate the BlastHitList for this query.
      BlastHitList *merged_bhl = Blast_HitListNew(1);

      // hsplist_count
      total_hsplist_count=0;//ceb Reset here for each query
      for ( int i=0; i < num_team_leaders; ++i )
	{
          BlastHitList *bhl = (hsp_results[i])->hitlist_array[BlastHSPResults_query_index];
          //dummy_int4 = bhl->hsplist_count;
	  //total_hsplist_count += dummy_int4;
	  //hsplist_counts[i] = dummy_int4;
	  total_hsplist_count += bhl->hsplist_count;
	  hsplist_counts[i] = bhl->hsplist_count;
	}

      merged_bhl->hsplist_count = total_hsplist_count;

      // hsplist_max
      total_hsplist_max=0;
      for ( int i=0; i < num_team_leaders; ++i )
	{
          BlastHitList *bhl = (hsp_results[i])->hitlist_array[BlastHSPResults_query_index];
	  if ( i == 0 )
	    {
	      total_hsplist_max = bhl->hsplist_max;
	    }
	}

      while ( total_hsplist_max < total_hsplist_count )
	{total_hsplist_max *= 2;}

      merged_bhl->hsplist_max = total_hsplist_max;

      // worst evalue
      worst_evalue_all = ((hsp_results[0])->hitlist_array[BlastHSPResults_query_index])->worst_evalue;
      for ( int i=1; i < num_team_leaders; ++i )
	{
          BlastHitList *bhl = (hsp_results[i])->hitlist_array[BlastHSPResults_query_index];
          dummy_double = bhl->worst_evalue;
	  if( dummy_double > worst_evalue_all )
	    {worst_evalue_all = dummy_double;}
	}
      merged_bhl->worst_evalue = worst_evalue_all;


      // best bit(low) score
      low_score_all = ((hsp_results[0])->hitlist_array[BlastHSPResults_query_index])->low_score;
      for ( int i=1; i < num_team_leaders; ++i )
	{
          BlastHitList *bhl = (hsp_results[i])->hitlist_array[BlastHSPResults_query_index];
	  dummy_int4 = bhl->low_score;
	  if( dummy_int4 < low_score_all )
	    {low_score_all = dummy_int4;}
	}
      merged_bhl->low_score = low_score_all;

      //heapified boolean
      all_heapified=0;
      for ( int i=0; i < num_team_leaders; ++i )
	{
          BlastHitList *bhl = (hsp_results[i])->hitlist_array[BlastHSPResults_query_index];
	  dummy_boolean = bhl->heapified;
	  if ( dummy_boolean == 1 )
	    {all_heapified = dummy_boolean;}
	}
      // if any are heapified assume the merge can be.
      merged_bhl->heapified = all_heapified;

      //hsplist_current
      dummy_int4 = ((hsp_results[0])->hitlist_array[BlastHSPResults_query_index])->hsplist_current;

      // take the first one or the max. add 100 until it surpasses hsplist_count!
      while ( dummy_int4 < total_hsplist_count )
	{dummy_int4 += 100;}

      total_hsplist_current = dummy_int4;
      merged_bhl->hsplist_current = total_hsplist_current;

      // Allocate the hspList to current long.
      merged_bhl->hsplist_array = (BlastHSPList**)calloc( merged_bhl->hsplist_current, sizeof(BlastHSPList*) );

      Int4 global_hsplist_index = 0;
      // Loop over each team leader and process the hit lists.
      for ( int i=0; i < num_team_leaders; ++i )
      {
        BlastHitList *bhl = (hsp_results[i])->hitlist_array[BlastHSPResults_query_index];
	for ( int BlastHitList_hsplist_index=0; BlastHitList_hsplist_index < hsplist_counts[i]; ++BlastHitList_hsplist_index )
	    {
	      // Allocate a BlastHSPList ptr
	      BlastHSPList *bhspl = bhl->hsplist_array[BlastHitList_hsplist_index];

#if USE_SHALLOW_COPY
	      BlastHSPList* merged_bhspl = (BlastHSPList*) calloc(1,sizeof(BlastHSPList) );

	      merged_bhspl->oid = bhspl->oid;
	      merged_bhspl->query_index = bhspl->query_index;
	      merged_bhspl->hspcnt = bhspl->hspcnt;
	      merged_bhspl->allocated = bhspl->allocated;
	      merged_bhspl->hsp_max = bhspl->hsp_max;
	      merged_bhspl->do_not_reallocate = bhspl->do_not_reallocate;
              merged_bhspl->best_evalue = bhspl->best_evalue;

	      // Allocate the array for High Scoring segment pairs.
	      merged_bhspl->hsp_array = (BlastHSP**) calloc( merged_bhspl->allocated, sizeof(BlastHSP*) );

	      // Loop over all HSPs between the query and the subject sequence.
	      for ( int BlastHSP_hsp_array_index=0; BlastHSP_hsp_array_index < bhspl->hspcnt; ++BlastHSP_hsp_array_index )
		{
		  BlastHSP* merged_bhsp = (BlastHSP*) calloc (1, sizeof(BlastHSP) );
		  BlastHSP* bhsp = bhspl->hsp_array[BlastHSP_hsp_array_index];

		  merged_bhsp->score = bhsp->score;
		  merged_bhsp->num_ident = bhsp->num_ident;
		  merged_bhsp->bit_score = bhsp->bit_score;
		  merged_bhsp->evalue = bhsp->evalue;
		  merged_bhsp->context = bhsp->context;
		  merged_bhsp->num = bhsp->num;
		  merged_bhsp->comp_adjustment_method = bhsp->comp_adjustment_method;
		  merged_bhsp->num_positives = bhsp->num_positives;
		  merged_bhsp->query.frame = bhsp->query.frame;
		  merged_bhsp->query.offset = bhsp->query.offset;
		  merged_bhsp->query.end = bhsp->query.end;
		  merged_bhsp->query.gapped_start = bhsp->query.gapped_start;
		  merged_bhsp->subject.frame = bhsp->subject.frame;
		  merged_bhsp->subject.offset = bhsp->subject.offset;
		  merged_bhsp->subject.end = bhsp->subject.end;
		  merged_bhsp->subject.gapped_start = bhsp->subject.gapped_start;

		  SPHIHspInfo *my_pat_info = bhsp->pat_info;
		  SPHIHspInfo *merged_pat_info = NULL;

                  if(my_pat_info!=NULL)
		  {
		    merged_pat_info = (SPHIHspInfo*) calloc( 1, sizeof(SPHIHspInfo) );
		    merged_pat_info->index = my_pat_info->index;
		    merged_pat_info->length = my_pat_info->length;
		  }
		  merged_bhsp->pat_info = merged_pat_info ;
		  
		  GapEditScript* my_gap_info = bhsp->gap_info;
		  GapEditScript *merged_gap_info = NULL;

                  if(my_gap_info != NULL)
		  {
		    merged_gap_info = (GapEditScript*) calloc(1, sizeof(GapEditScript) );
		    merged_gap_info->size = my_gap_info->size;
		    merged_gap_info->op_type = (EGapAlignOpType*) calloc( merged_gap_info->size , sizeof(EGapAlignOpType) );
		    memcpy(&(merged_gap_info->op_type[0]) , &(my_gap_info->op_type[0]) , my_gap_info->size * sizeof(EGapAlignOpType) );
		    merged_gap_info->num = (Int4*) calloc( merged_gap_info->size , sizeof(Int4) );
		    memcpy( &(merged_gap_info->num[0]) , &(my_gap_info->num[0]) , my_gap_info->size * sizeof(Int4) );
		  }
		  merged_bhsp->gap_info = merged_gap_info;

		  // Now attach this BlastHSP to the parent BlastHSPList.
		  merged_bhspl->hsp_array[ BlastHSP_hsp_array_index ] = merged_bhsp ;

		} // End of loop over the BlastHSP_hsp_array_index -> the high scoring pair data


#else
	      // Now attach this BlastHSP to the parent BlastHSPList.
	      BlastHSPList* merged_bhspl = BlastHSPListDuplicate(bhspl);
#endif


	      // Now attach this BlastHSPList to the parent BlastHitList
	      merged_bhl->hsplist_array[ global_hsplist_index ] = merged_bhspl ;
	      ++global_hsplist_index;  // increment the global index
	      
	    }  // End loop over BlastHitList_hsplist_index

	} // End Loop over team leaders.


      // Now we can add this BlastHitList to the parent BlastHSPResults
      merged_hsp_results->hitlist_array[ BlastHSPResults_query_index ] = merged_bhl ;
  }

  // Now we can sort the HSPResults.
  Blast_HSPResultsSortByEvalue( merged_hsp_results ) ;

  return merged_hsp_results;
}



BlastHSPList* BlastHSPListDuplicate(const BlastHSPList *bhspl)
{
  // Allocate a BlastHSPList ptr
#if 0
  BlastHSPList* new_bhspl = BlastHSPListDup(bhspl);

#else  
  BlastHSPList* new_bhspl = Blast_HSPListNew(bhspl->hsp_max);
  new_bhspl->oid = bhspl->oid;
  new_bhspl->query_index = bhspl->query_index;
  new_bhspl->hspcnt = bhspl->hspcnt;
  new_bhspl->allocated = bhspl->allocated;
  new_bhspl->hsp_max = bhspl->hsp_max;
  new_bhspl->do_not_reallocate = bhspl->do_not_reallocate;
  new_bhspl->best_evalue = bhspl->best_evalue;
 
  // Allocate the array for High Scoring segment pairs.
  //new_bhspl->hsp_array = (BlastHSP**) calloc( new_bhspl->allocated, sizeof(BlastHSP*) );
  //new_bhspl->hsp_array = Blast_HSPListNew(new_bhspl->hsp_max);

  new_bhspl->hsp_array = (BlastHSP**)realloc( new_bhspl->hsp_array,new_bhspl->allocated*sizeof(BlastHSP*));

  // Loop over all HSPs between the query and the subject sequence.
  for ( int BlastHSP_hsp_array_index=0; BlastHSP_hsp_array_index < bhspl->hspcnt; ++BlastHSP_hsp_array_index )
    {
      //BlastHSP* new_bhsp = (BlastHSP*) calloc (1, sizeof(BlastHSP) );
      //BlastHSP* new_bhsp = Blast_HSPNew();
      BlastHSP* bhsp = bhspl->hsp_array[BlastHSP_hsp_array_index];
#if 1
      //works
      BlastHSP* new_bhsp = NULL;
      /* Do not pass the edit script, because we don't want to tranfer 
	 ownership. */
      Blast_HSPInit(bhsp->query.offset, bhsp->query.end, bhsp->subject.offset, 
		    bhsp->subject.end, bhsp->query.gapped_start, 
		    bhsp->subject.gapped_start, bhsp->context, 
		    bhsp->query.frame, bhsp->subject.frame, bhsp->score, 
		    NULL, &new_bhsp);
      new_bhsp->evalue = bhsp->evalue;
      new_bhsp->num = bhsp->num;
      new_bhsp->num_ident = bhsp->num_ident;
      new_bhsp->bit_score = bhsp->bit_score;
      new_bhsp->comp_adjustment_method = bhsp->comp_adjustment_method;
      if (bhsp->gap_info) {
	new_bhsp->gap_info = GapEditScriptDup(bhsp->gap_info);
      }
      
      if (bhsp->pat_info) {
	/* Copy this HSP's pattern data. */
	new_bhsp->pat_info = 
	  (SPHIHspInfo*) BlastMemDup(bhsp->pat_info, sizeof(SPHIHspInfo));
      }
      
#else
      BlastHSP* new_bhsp = Blast_HSPNew();
      
      new_bhsp->score = bhsp->score;
      new_bhsp->num_ident = bhsp->num_ident;
      new_bhsp->bit_score = bhsp->bit_score;
      new_bhsp->evalue = bhsp->evalue;
      new_bhsp->context = bhsp->context;
      new_bhsp->num = bhsp->num;
      new_bhsp->comp_adjustment_method = bhsp->comp_adjustment_method;
      new_bhsp->num_positives = bhsp->num_positives;
      new_bhsp->query.frame = bhsp->query.frame;
      new_bhsp->query.offset = bhsp->query.offset;
      new_bhsp->query.end = bhsp->query.end;
      new_bhsp->query.gapped_start = bhsp->query.gapped_start;
      new_bhsp->subject.frame = bhsp->subject.frame;
      new_bhsp->subject.offset = bhsp->subject.offset;
      new_bhsp->subject.end = bhsp->subject.end;
      new_bhsp->subject.gapped_start = bhsp->subject.gapped_start;
      
      SPHIHspInfo *my_pat_info = bhsp->pat_info;
      SPHIHspInfo *new_pat_info = NULL;
      
      if(my_pat_info!=NULL)
	{
	  new_pat_info = (SPHIHspInfo*) calloc( 1, sizeof(SPHIHspInfo) );
	  new_pat_info->index = my_pat_info->index;
	  new_pat_info->length = my_pat_info->length;
	}
      new_bhsp->pat_info = new_pat_info ;
      
      GapEditScript* my_gap_info = bhsp->gap_info;
      GapEditScript *new_gap_info = NULL;
      
      if(my_gap_info != NULL)
	{
	  new_gap_info = (GapEditScript*) calloc(1, sizeof(GapEditScript) );
	  new_gap_info->size = my_gap_info->size;
	  new_gap_info->op_type = (EGapAlignOpType*) calloc( new_gap_info->size , sizeof(EGapAlignOpType) );
	  memcpy(&(new_gap_info->op_type[0]) , &(my_gap_info->op_type[0]) , my_gap_info->size * sizeof(EGapAlignOpType) );
	  new_gap_info->num = (Int4*) calloc( new_gap_info->size , sizeof(Int4) );
	  memcpy( &(new_gap_info->num[0]) , &(my_gap_info->num[0]) , my_gap_info->size * sizeof(Int4) );
	}
      new_bhsp->gap_info = new_gap_info;
#endif

      // Now attach this BlastHSP to the parent BlastHSPList.
      new_bhspl->hsp_array[ BlastHSP_hsp_array_index ] = new_bhsp ;
      
    } // End of loop over the BlastHSP_hsp_array_index -> the high scoring pair data
#endif  

  return new_bhspl;  
}


//======================================================================




BlastHSPResults* gbhrpfncDuplicateBlastHSPResults ( BlastHSPResults* HSPResultsArray)
{
  int tid = 0;
  int num_team_leaders=1; 
  Int4 num_queries = 0;
  Int4 dummy_int4 = 0;  // dummy variable.
  double dummy_double = 0.;
  Boolean dummy_boolean;
  Boolean all_heapified = 0;

  Int4 total_hsplist_current = 0;
  Int4 total_hsplist_count = 0;
  Int4 total_hsplist_max=0;
  Int4 total_hsparray_allocated = 0;

  double worst_evalue_all=0.;
  Int4 low_score_all=0;

  // Need a buffer pointer for each different buffer.
  Int4 *hsplist_counts = NULL;
  BlastHSPResults **hsp_results =  new BlastHSPResults*[num_team_leaders];    // use this to point to a value in HSPResultsArray
  BlastHSPResults *merged_hsp_results = NULL;    // The guy we're returning.

  hsplist_counts = (Int4*)malloc(num_team_leaders * sizeof(char*));

  for ( int i=0; i < num_team_leaders; ++i )
  {
    hsp_results[i] = HSPResultsArray;
  }

  // Get the number of queries -- its the same across all threads in a thread group.
  num_queries = (hsp_results[0])->num_queries;
  //std::cout<<"(hsp_results[0])->num_queries= "<<num_queries <<std::endl;std::cout.flush();
  merged_hsp_results = Blast_HSPResultsNew(num_queries);

  // We have to move across the buffer in a very particular way (reverse of how its packed.)
  // See traceback_stage.cpp for details on how its packed.

  for ( int BlastHSPResults_query_index = 0; BlastHSPResults_query_index < merged_hsp_results->num_queries ; ++BlastHSPResults_query_index )
    {    
      // Allocate the BlastHitList for this query.
      BlastHitList *merged_bhl = Blast_HitListNew(1);

      // hsplist_count
      total_hsplist_count=0;//ceb Reset here for each query
      for ( int i=0; i < num_team_leaders; ++i )
	{
          BlastHitList *bhl = (hsp_results[i])->hitlist_array[BlastHSPResults_query_index];
          dummy_int4 = bhl->hsplist_count;
	  total_hsplist_count += dummy_int4;
	  hsplist_counts[i] = dummy_int4;
	}

      merged_bhl->hsplist_count = total_hsplist_count;

      // hsplist_max
      total_hsplist_max=0;
      for ( int i=0; i < num_team_leaders; ++i )
	{
          BlastHitList *bhl = (hsp_results[i])->hitlist_array[BlastHSPResults_query_index];
	  if ( i == 0 )
	    {
	      total_hsplist_max = bhl->hsplist_max;
	    }
	}

      while ( total_hsplist_max < total_hsplist_count )
	{total_hsplist_max *= 2;}

      merged_bhl->hsplist_max = total_hsplist_max;

      // worst evalue
      worst_evalue_all = ((hsp_results[0])->hitlist_array[BlastHSPResults_query_index])->worst_evalue;
      for ( int i=1; i < num_team_leaders; ++i )
	{
          BlastHitList *bhl = (hsp_results[i])->hitlist_array[BlastHSPResults_query_index];
          dummy_double = bhl->worst_evalue;
	  if( dummy_double > worst_evalue_all )
	    {worst_evalue_all = dummy_double;}
	}
      merged_bhl->worst_evalue = worst_evalue_all;


      // best bit(low) score
      low_score_all = ((hsp_results[0])->hitlist_array[BlastHSPResults_query_index])->low_score;
      for ( int i=1; i < num_team_leaders; ++i )
	{
          BlastHitList *bhl = (hsp_results[i])->hitlist_array[BlastHSPResults_query_index];
	  dummy_int4 = bhl->low_score;
	  if( dummy_int4 < low_score_all )
	    {low_score_all = dummy_int4;}
	}
      merged_bhl->low_score = low_score_all;

      //heapified boolean
      all_heapified=0;
      for ( int i=0; i < num_team_leaders; ++i )
	{
          BlastHitList *bhl = (hsp_results[i])->hitlist_array[BlastHSPResults_query_index];
	  dummy_boolean = bhl->heapified;
	  if ( dummy_boolean == 1 )
	    {all_heapified = dummy_boolean;}
	}
      // if any are heapified assume the merge can be.
      merged_bhl->heapified = all_heapified;

      //hsplist_current
      dummy_int4 = ((hsp_results[0])->hitlist_array[BlastHSPResults_query_index])->hsplist_current;

      // take the first one or the max. add 100 until it surpasses hsplist_count!
      while ( dummy_int4 < total_hsplist_count )
	{dummy_int4 += 100;}

      total_hsplist_current = dummy_int4;
      merged_bhl->hsplist_current = total_hsplist_current;

      // Allocate the hspList to current long.
      merged_bhl->hsplist_array = (BlastHSPList**)calloc( merged_bhl->hsplist_current, sizeof(BlastHSPList*) );

      Int4 global_hsplist_index = 0;
      // Loop over each team leader and process the hit lists.
      for ( int i=0; i < num_team_leaders; ++i )
      {
        BlastHitList *bhl = (hsp_results[i])->hitlist_array[BlastHSPResults_query_index];
	for ( int BlastHitList_hsplist_index=0; BlastHitList_hsplist_index < hsplist_counts[i]; ++BlastHitList_hsplist_index )
	    {
	      // Allocate a BlastHSPList ptr
	      BlastHSPList* merged_bhspl = (BlastHSPList*) calloc(1,sizeof(BlastHSPList) );
	      BlastHSPList *bhspl = bhl->hsplist_array[BlastHitList_hsplist_index];

	      merged_bhspl->oid = bhspl->oid;
	      merged_bhspl->query_index = bhspl->query_index;
	      merged_bhspl->hspcnt = bhspl->hspcnt;
	      merged_bhspl->allocated = bhspl->allocated;
	      merged_bhspl->hsp_max = bhspl->hsp_max;
	      merged_bhspl->do_not_reallocate = bhspl->do_not_reallocate;
              merged_bhspl->best_evalue = bhspl->best_evalue;

	      // Allocate the array for High Scoring segment pairs.
	      merged_bhspl->hsp_array = (BlastHSP**) calloc( merged_bhspl->allocated, sizeof(BlastHSP*) );

	      // Loop over all HSPs between the query and the subject sequence.
	      for ( int BlastHSP_hsp_array_index=0; BlastHSP_hsp_array_index < bhspl->hspcnt; ++BlastHSP_hsp_array_index )
		{
		  BlastHSP* merged_bhsp = (BlastHSP*) calloc (1, sizeof(BlastHSP) );
		  BlastHSP* bhsp = bhspl->hsp_array[BlastHSP_hsp_array_index];

		  merged_bhsp->score = bhsp->score;
		  merged_bhsp->num_ident = bhsp->num_ident;
		  merged_bhsp->bit_score = bhsp->bit_score;
		  merged_bhsp->evalue = bhsp->evalue;
		  merged_bhsp->context = bhsp->context;
		  merged_bhsp->num = bhsp->num;
		  merged_bhsp->comp_adjustment_method = bhsp->comp_adjustment_method;
		  merged_bhsp->num_positives = bhsp->num_positives;
		  merged_bhsp->query.frame = bhsp->query.frame;
		  merged_bhsp->query.offset = bhsp->query.offset;
		  merged_bhsp->query.end = bhsp->query.end;
		  merged_bhsp->query.gapped_start = bhsp->query.gapped_start;
		  merged_bhsp->subject.frame = bhsp->subject.frame;
		  merged_bhsp->subject.offset = bhsp->subject.offset;
		  merged_bhsp->subject.end = bhsp->subject.end;
		  merged_bhsp->subject.gapped_start = bhsp->subject.gapped_start;

		  SPHIHspInfo *my_pat_info = bhsp->pat_info;
		  SPHIHspInfo *merged_pat_info = NULL;

                  if(my_pat_info!=NULL)
		  {
		    merged_pat_info = (SPHIHspInfo*) calloc( 1, sizeof(SPHIHspInfo) );
		    merged_pat_info->index = my_pat_info->index;
		    merged_pat_info->length = my_pat_info->length;
		  }
		  merged_bhsp->pat_info = merged_pat_info ;
		  
		  GapEditScript* my_gap_info = bhsp->gap_info;
		  GapEditScript *merged_gap_info = NULL;

                  if(my_gap_info != NULL)
		  {
		    merged_gap_info = (GapEditScript*) calloc(1, sizeof(GapEditScript) );
		    merged_gap_info->size = my_gap_info->size;
		    merged_gap_info->op_type = (EGapAlignOpType*) calloc( merged_gap_info->size , sizeof(EGapAlignOpType) );
		    memcpy(&(merged_gap_info->op_type[0]) , &(my_gap_info->op_type[0]) , my_gap_info->size * sizeof(EGapAlignOpType) );
		    merged_gap_info->num = (Int4*) calloc( merged_gap_info->size , sizeof(Int4) );
		    memcpy( &(merged_gap_info->num[0]) , &(my_gap_info->num[0]) , my_gap_info->size * sizeof(Int4) );
		  }
		  merged_bhsp->gap_info = merged_gap_info;

		  // Now attach this BlastHSP to the parent BlastHSPList.
		  merged_bhspl->hsp_array[ BlastHSP_hsp_array_index ] = merged_bhsp ;

		} // End of loop over the BlastHSP_hsp_array_index -> the high scoring pair data

	      // Now attach this BlastHSPList to the parent BlastHitList
	      merged_bhl->hsplist_array[ global_hsplist_index ] = merged_bhspl ;
	      ++global_hsplist_index;  // increment the global index
	      
	    }  // End loop over BlastHitList_hsplist_index

	} // End Loop over team leaders.


      // Now we can add this BlastHitList to the parent BlastHSPResults
      merged_hsp_results->hitlist_array[ BlastHSPResults_query_index ] = merged_bhl ;
  }

  // Now we can sort the HSPResults.
  Blast_HSPResultsSortByEvalue( merged_hsp_results ) ;

  return merged_hsp_results;
}





BlastScoreBlk *gbsbpfncDuplicateBlastScoreBlk(BlastScoreBlk *ScoreBlkGroup)
{

#if 0
  for (index=0; index<sbp->number_of_contexts; index++) {
    if (sbp->sfp)
      sbp->sfp[index] = Blast_ScoreFreqFree(sbp->sfp[index]);
    if (sbp->kbp_std)
      sbp->kbp_std[index] = Blast_KarlinBlkFree(sbp->kbp_std[index]);
    if (sbp->kbp_gap_std)
      sbp->kbp_gap_std[index] = Blast_KarlinBlkFree(sbp->kbp_gap_std[index]);
    if (sbp->kbp_psi)
      sbp->kbp_psi[index] = Blast_KarlinBlkFree(sbp->kbp_psi[index]);
    if (sbp->kbp_gap_psi)
      sbp->kbp_gap_psi[index] = Blast_KarlinBlkFree(sbp->kbp_gap_psi[index]);
  }
  if (sbp->kbp_ideal)
    sbp->kbp_ideal = Blast_KarlinBlkFree(sbp->kbp_ideal);
  if (sbp->gbp) 
    sbp->gbp = s_BlastGumbelBlkFree(sbp->gbp);
  sfree(sbp->sfp);
  sfree(sbp->kbp_std);
  sfree(sbp->kbp_psi);
  sfree(sbp->kbp_gap_std);
  sfree(sbp->kbp_gap_psi);
  sbp->matrix = SBlastScoreMatrixFree(sbp->matrix);
  sbp->comments = ListNodeFreeData(sbp->comments);
  sfree(sbp->name);
  sbp->psi_matrix = SPsiBlastScoreMatrixFree(sbp->psi_matrix);
  sfree(sbp->ambiguous_res);
  sfree(sbp);
#endif


  if(ScoreBlkGroup==NULL){return NULL;}

  BlastScoreBlk* lbsbpBlastScoreBlk = (BlastScoreBlk*)calloc(1,sizeof(BlastScoreBlk));

  lbsbpBlastScoreBlk->protein_alphabet=ScoreBlkGroup->protein_alphabet;

  lbsbpBlastScoreBlk->alphabet_code=ScoreBlkGroup->alphabet_code;

  lbsbpBlastScoreBlk->alphabet_size=ScoreBlkGroup->alphabet_size;

  lbsbpBlastScoreBlk->alphabet_start=ScoreBlkGroup->alphabet_start;

  //int str_length = strlen( ScoreBlkGroup->name );
  //lbsbpBlastScoreBlk->name = malloc( sizeof(char)*(str_length+1) );
  //strncpy ( lbsbpBlastScoreBlk->name, ScoreBlkGroup->name, (size_t) str_length );
  lbsbpBlastScoreBlk->name = ScoreBlkGroup->name;

  lbsbpBlastScoreBlk->comments=NULL;

  ListNode *llndpCurrentNode=NULL;

  ListNode *llndpPreviousNode=NULL;


 {

    if (ScoreBlkGroup->comments) {

      if (lbsbpBlastScoreBlk->comments) {

        llndpPreviousNode=lbsbpBlastScoreBlk->comments;

        llndpCurrentNode=lbsbpBlastScoreBlk->comments->next;

         while (llndpCurrentNode) {

           llndpPreviousNode=llndpCurrentNode;

           llndpCurrentNode=llndpCurrentNode->next;

         }

         llndpPreviousNode->next=ScoreBlkGroup->comments;

       } else {

         lbsbpBlastScoreBlk->comments=ScoreBlkGroup->comments;

      }

    }

  }


  lbsbpBlastScoreBlk->matrix = NULL;

 {

    if (ScoreBlkGroup->matrix) {

      lbsbpBlastScoreBlk->matrix=ScoreBlkGroup->matrix;
      //ceb make deep copy

    }

  }


  lbsbpBlastScoreBlk->psi_matrix = NULL;

  {
    if (ScoreBlkGroup->psi_matrix) {
      lbsbpBlastScoreBlk->psi_matrix=ScoreBlkGroup->psi_matrix;
    }
  }


  lbsbpBlastScoreBlk->matrix_only_scoring=ScoreBlkGroup->matrix_only_scoring;
  lbsbpBlastScoreBlk->complexity_adjusted_scoring=ScoreBlkGroup->complexity_adjusted_scoring;
  lbsbpBlastScoreBlk->loscore=ScoreBlkGroup->loscore;
  lbsbpBlastScoreBlk->hiscore=ScoreBlkGroup->hiscore;
  lbsbpBlastScoreBlk->penalty=ScoreBlkGroup->penalty;
  lbsbpBlastScoreBlk->reward=ScoreBlkGroup->reward;
  lbsbpBlastScoreBlk->scale_factor=ScoreBlkGroup->scale_factor;
  lbsbpBlastScoreBlk->read_in_matrix=ScoreBlkGroup->read_in_matrix;
  


  int lintTotalNumberOfContexts=0;

  {
    lintTotalNumberOfContexts+=ScoreBlkGroup->number_of_contexts;
    //lintTotalNumberOfContexts+=ScoreBlkGroup->number_of_contexts;
  }

  lbsbpBlastScoreBlk->number_of_contexts=lintTotalNumberOfContexts;

  lbsbpBlastScoreBlk->sfp=NULL;

  if (lbsbpBlastScoreBlk->number_of_contexts>0) { 

    if (ScoreBlkGroup->sfp) {

      lbsbpBlastScoreBlk->sfp = (Blast_ScoreFreq**) calloc( lbsbpBlastScoreBlk->number_of_contexts, sizeof(Blast_ScoreFreq*) );

      int lintGlobalContextNumber=-1;

      {

        for (Int4 j=0; j < ScoreBlkGroup->number_of_contexts; ++j) {

          lintGlobalContextNumber++;

          lbsbpBlastScoreBlk->sfp[lintGlobalContextNumber]=ScoreBlkGroup->sfp[j];

        }

      }

    }

  }



  lbsbpBlastScoreBlk->kbp=NULL;

  if (lintTotalNumberOfContexts>0) {
    if (ScoreBlkGroup->kbp) {
      lbsbpBlastScoreBlk->kbp = (Blast_KarlinBlk**) calloc( lbsbpBlastScoreBlk->number_of_contexts, sizeof(Blast_KarlinBlk*) );
      int lintGlobalContextNumber=-1;

       {
        for (Int4 j=0; j < ScoreBlkGroup->number_of_contexts; ++j) {
          lintGlobalContextNumber++;
          lbsbpBlastScoreBlk->kbp[lintGlobalContextNumber]=ScoreBlkGroup->kbp[j];
        }
      }
    }
  }



  lbsbpBlastScoreBlk->kbp_gap=NULL;

  if (lintTotalNumberOfContexts>0) {
    if (ScoreBlkGroup->kbp_gap) {
      lbsbpBlastScoreBlk->kbp_gap = (Blast_KarlinBlk**) calloc( lbsbpBlastScoreBlk->number_of_contexts, sizeof(Blast_KarlinBlk*) );
      int lintGlobalContextNumber=-1;

      {
        for (Int4 j=0; j < ScoreBlkGroup->number_of_contexts; ++j) {
          lintGlobalContextNumber++;
          lbsbpBlastScoreBlk->kbp_gap[lintGlobalContextNumber]=ScoreBlkGroup->kbp_gap[j];
        }
      }
    }
  }


  lbsbpBlastScoreBlk->gbp = NULL;

  {
    if (ScoreBlkGroup->gbp) {
      lbsbpBlastScoreBlk->gbp=ScoreBlkGroup->gbp;
    }
  }


  lbsbpBlastScoreBlk->kbp_std=NULL;

  if (lintTotalNumberOfContexts>0) {
    if (ScoreBlkGroup->kbp_std) {
      lbsbpBlastScoreBlk->kbp_std = (Blast_KarlinBlk**) calloc( lbsbpBlastScoreBlk->number_of_contexts, sizeof(Blast_KarlinBlk*) );
      int lintGlobalContextNumber=-1;

      {
        for (Int4 j=0; j < ScoreBlkGroup->number_of_contexts; ++j) {
          lintGlobalContextNumber++;
          lbsbpBlastScoreBlk->kbp_std[lintGlobalContextNumber]=ScoreBlkGroup->kbp_std[j];
        }
      }
    }
  }

  lbsbpBlastScoreBlk->kbp_psi=NULL;

  if (lintTotalNumberOfContexts>0) {
    if (ScoreBlkGroup->kbp_psi) {
      lbsbpBlastScoreBlk->kbp_psi = (Blast_KarlinBlk**) calloc( lbsbpBlastScoreBlk->number_of_contexts, sizeof(Blast_KarlinBlk*) );
      int lintGlobalContextNumber=-1;

       {
        for (Int4 j=0; j < ScoreBlkGroup->number_of_contexts; ++j ) {
          lintGlobalContextNumber++;
          lbsbpBlastScoreBlk->kbp_psi[lintGlobalContextNumber] = ScoreBlkGroup->kbp_psi[j];
        }
      }
    }
  }

  lbsbpBlastScoreBlk->kbp_gap_std = NULL;

  if (lintTotalNumberOfContexts>0) {
    if (ScoreBlkGroup->kbp_gap_std) {
      lbsbpBlastScoreBlk->kbp_gap_std = (Blast_KarlinBlk**) calloc( lbsbpBlastScoreBlk->number_of_contexts, sizeof(Blast_KarlinBlk*) );
      int lintGlobalContextNumber = -1;

      {
        for (Int4 j=0; j < ScoreBlkGroup->number_of_contexts; ++j ) 
        {
          lintGlobalContextNumber++;
          lbsbpBlastScoreBlk->kbp_gap_std[lintGlobalContextNumber] = ScoreBlkGroup->kbp_gap_std[j];
        }
      }
    }
  }

  lbsbpBlastScoreBlk->kbp_gap_psi=NULL;

  if (lintTotalNumberOfContexts>0) 
  {
    if (ScoreBlkGroup->kbp_gap_psi) 
    {
      lbsbpBlastScoreBlk->kbp_gap_psi = (Blast_KarlinBlk**) calloc( lbsbpBlastScoreBlk->number_of_contexts, sizeof(Blast_KarlinBlk*) );
      int lintGlobalContextNumber = -1;


      {
        for (Int4 j=0; j < ScoreBlkGroup->number_of_contexts; ++j ) 
        {
          lintGlobalContextNumber++;
          lbsbpBlastScoreBlk->kbp_gap_psi[lintGlobalContextNumber] = ScoreBlkGroup->kbp_gap_psi[j];
        }
      }
    }
  }

  lbsbpBlastScoreBlk->kbp_ideal = NULL;

  {
    if (ScoreBlkGroup->kbp_ideal) 
    {
      lbsbpBlastScoreBlk->kbp_ideal = ScoreBlkGroup->kbp_ideal;
    }
  }

  lbsbpBlastScoreBlk->ambig_size = ScoreBlkGroup->ambig_size;
  lbsbpBlastScoreBlk->ambig_occupy = ScoreBlkGroup->ambig_occupy;
  lbsbpBlastScoreBlk->ambiguous_res = ScoreBlkGroup->ambiguous_res;
  lbsbpBlastScoreBlk->round_down = ScoreBlkGroup->round_down;

  return lbsbpBlastScoreBlk;
}




// Functions to combine objects across thread groups

BlastScoreBlk *gbsbpfncCombineThreadGroupBlastScoreBlk(int aintNumberOfThreadGroups, BlastScoreBlk **ScoreBlkGroup) {

  BlastScoreBlk* lbsbpBlastScoreBlk = (BlastScoreBlk*)calloc(1,sizeof(BlastScoreBlk));

  lbsbpBlastScoreBlk->protein_alphabet=ScoreBlkGroup[0]->protein_alphabet;

  lbsbpBlastScoreBlk->alphabet_code=ScoreBlkGroup[0]->alphabet_code;

  lbsbpBlastScoreBlk->alphabet_size=ScoreBlkGroup[0]->alphabet_size;

  lbsbpBlastScoreBlk->alphabet_start=ScoreBlkGroup[0]->alphabet_start;

  lbsbpBlastScoreBlk->name=ScoreBlkGroup[0]->name;

  lbsbpBlastScoreBlk->comments=NULL;

  ListNode *llndpCurrentNode=NULL;

  ListNode *llndpPreviousNode=NULL;


  for (int i=0; i<aintNumberOfThreadGroups; i++) {

    if (ScoreBlkGroup[i]->comments) {

      if (lbsbpBlastScoreBlk->comments) {

        llndpPreviousNode=lbsbpBlastScoreBlk->comments;

        llndpCurrentNode=lbsbpBlastScoreBlk->comments->next;

         while (llndpCurrentNode) {

           llndpPreviousNode=llndpCurrentNode;

           llndpCurrentNode=llndpCurrentNode->next;

         }

         llndpPreviousNode->next=ScoreBlkGroup[i]->comments;

       } else {

         lbsbpBlastScoreBlk->comments=ScoreBlkGroup[i]->comments;

      }

    }

  }


  lbsbpBlastScoreBlk->matrix = NULL;

  for (int i=0; i<aintNumberOfThreadGroups; i++) {

    if (ScoreBlkGroup[i]->matrix) {

      lbsbpBlastScoreBlk->matrix=ScoreBlkGroup[i]->matrix;

    }

  }


  lbsbpBlastScoreBlk->psi_matrix = NULL;

  for (int i=0; i<aintNumberOfThreadGroups; i++) {
    if (ScoreBlkGroup[i]->psi_matrix) {
      lbsbpBlastScoreBlk->psi_matrix=ScoreBlkGroup[i]->psi_matrix;
    }
  }


  lbsbpBlastScoreBlk->matrix_only_scoring=ScoreBlkGroup[0]->matrix_only_scoring;
  lbsbpBlastScoreBlk->complexity_adjusted_scoring=ScoreBlkGroup[0]->complexity_adjusted_scoring;
  lbsbpBlastScoreBlk->loscore=ScoreBlkGroup[0]->loscore;
  lbsbpBlastScoreBlk->hiscore=ScoreBlkGroup[0]->hiscore;
  lbsbpBlastScoreBlk->penalty=ScoreBlkGroup[0]->penalty;
  lbsbpBlastScoreBlk->reward=ScoreBlkGroup[0]->reward;
  lbsbpBlastScoreBlk->scale_factor=ScoreBlkGroup[0]->scale_factor;
  lbsbpBlastScoreBlk->read_in_matrix=ScoreBlkGroup[0]->read_in_matrix;
  


  int lintTotalNumberOfContexts=0;

  for (int i=0; i<aintNumberOfThreadGroups; i++) {
    lintTotalNumberOfContexts+=ScoreBlkGroup[i]->number_of_contexts;
    //lintTotalNumberOfContexts+=ScoreBlkGroup[0]->number_of_contexts;
  }

  lbsbpBlastScoreBlk->number_of_contexts=lintTotalNumberOfContexts;

  lbsbpBlastScoreBlk->sfp=NULL;

  if (lbsbpBlastScoreBlk->number_of_contexts>0) { 

    if (ScoreBlkGroup[0]->sfp) {

      lbsbpBlastScoreBlk->sfp = (Blast_ScoreFreq**) calloc( lbsbpBlastScoreBlk->number_of_contexts, sizeof(Blast_ScoreFreq*) );

      int lintGlobalContextNumber=-1;

      for (int i=0; i<aintNumberOfThreadGroups; i++) {

        for (Int4 j=0; j < ScoreBlkGroup[i]->number_of_contexts; ++j) {

          lintGlobalContextNumber++;

          lbsbpBlastScoreBlk->sfp[lintGlobalContextNumber]=ScoreBlkGroup[i]->sfp[j];

        }

      }

    }

  }



  lbsbpBlastScoreBlk->kbp=NULL;

  if (lintTotalNumberOfContexts>0) {
    if (ScoreBlkGroup[0]->kbp) {
      lbsbpBlastScoreBlk->kbp = (Blast_KarlinBlk**) calloc( lbsbpBlastScoreBlk->number_of_contexts, sizeof(Blast_KarlinBlk*) );
      int lintGlobalContextNumber=-1;

      for (int i=0; i<aintNumberOfThreadGroups; i++) {
        for (Int4 j=0; j < ScoreBlkGroup[i]->number_of_contexts; ++j) {
          lintGlobalContextNumber++;
          lbsbpBlastScoreBlk->kbp[lintGlobalContextNumber]=ScoreBlkGroup[i]->kbp[j];
        }
      }
    }
  }



  lbsbpBlastScoreBlk->kbp_gap=NULL;

  if (lintTotalNumberOfContexts>0) {
    if (ScoreBlkGroup[0]->kbp_gap) {
      lbsbpBlastScoreBlk->kbp_gap = (Blast_KarlinBlk**) calloc( lbsbpBlastScoreBlk->number_of_contexts, sizeof(Blast_KarlinBlk*) );
      int lintGlobalContextNumber=-1;

      for (int i=0; i<aintNumberOfThreadGroups; i++) {
        for (Int4 j=0; j < ScoreBlkGroup[i]->number_of_contexts; ++j) {
          lintGlobalContextNumber++;
          lbsbpBlastScoreBlk->kbp_gap[lintGlobalContextNumber]=ScoreBlkGroup[i]->kbp_gap[j];
        }
      }
    }
  }


  lbsbpBlastScoreBlk->gbp = NULL;

  for (int i=0; i<aintNumberOfThreadGroups; i++) {
    if (ScoreBlkGroup[i]->gbp) {
      lbsbpBlastScoreBlk->gbp=ScoreBlkGroup[i]->gbp;
    }
  }



  lbsbpBlastScoreBlk->kbp_std=NULL;

  if (lintTotalNumberOfContexts>0) {
    if (ScoreBlkGroup[0]->kbp_std) {
      lbsbpBlastScoreBlk->kbp_std = (Blast_KarlinBlk**) calloc( lbsbpBlastScoreBlk->number_of_contexts, sizeof(Blast_KarlinBlk*) );
      int lintGlobalContextNumber=-1;

      for (int i=0; i<aintNumberOfThreadGroups; i++) {
        for (Int4 j=0; j < ScoreBlkGroup[i]->number_of_contexts; ++j) {
          lintGlobalContextNumber++;
          lbsbpBlastScoreBlk->kbp_std[lintGlobalContextNumber]=ScoreBlkGroup[i]->kbp_std[j];
        }
      }
    }
  }



  lbsbpBlastScoreBlk->kbp_psi=NULL;

  if (lintTotalNumberOfContexts>0) {
    if (ScoreBlkGroup[0]->kbp_psi) {
      lbsbpBlastScoreBlk->kbp_psi = (Blast_KarlinBlk**) calloc( lbsbpBlastScoreBlk->number_of_contexts, sizeof(Blast_KarlinBlk*) );
      int lintGlobalContextNumber=-1;

      for (int i=0; i<aintNumberOfThreadGroups; i++) {
        for (Int4 j=0; j < ScoreBlkGroup[i]->number_of_contexts; ++j ) {
          lintGlobalContextNumber++;
          lbsbpBlastScoreBlk->kbp_psi[lintGlobalContextNumber] = ScoreBlkGroup[i]->kbp_psi[j];
        }
      }
    }
  }


  lbsbpBlastScoreBlk->kbp_gap_std = NULL;

  if (lintTotalNumberOfContexts>0) {
    if (ScoreBlkGroup[0]->kbp_gap_std) {
      lbsbpBlastScoreBlk->kbp_gap_std = (Blast_KarlinBlk**) calloc( lbsbpBlastScoreBlk->number_of_contexts, sizeof(Blast_KarlinBlk*) );
      int lintGlobalContextNumber = -1;

      for (int i=0; i<aintNumberOfThreadGroups; i++) 
      {
        for (Int4 j=0; j < ScoreBlkGroup[i]->number_of_contexts; ++j ) 
        {
          lintGlobalContextNumber++;
          lbsbpBlastScoreBlk->kbp_gap_std[lintGlobalContextNumber] = ScoreBlkGroup[i]->kbp_gap_std[j];
        }
      }
    }
  }


  lbsbpBlastScoreBlk->kbp_gap_psi=NULL;

  if (lintTotalNumberOfContexts>0) 
  {
    if (ScoreBlkGroup[0]->kbp_gap_psi) 
    {
      lbsbpBlastScoreBlk->kbp_gap_psi = (Blast_KarlinBlk**) calloc( lbsbpBlastScoreBlk->number_of_contexts, sizeof(Blast_KarlinBlk*) );
      int lintGlobalContextNumber = -1;

      for (int i=0; i < aintNumberOfThreadGroups; i++) 
      {
        for (Int4 j=0; j < ScoreBlkGroup[i]->number_of_contexts; ++j ) 
        {
          lintGlobalContextNumber++;
          lbsbpBlastScoreBlk->kbp_gap_psi[lintGlobalContextNumber] = ScoreBlkGroup[i]->kbp_gap_psi[j];
        }
      }
    }
  }


  lbsbpBlastScoreBlk->kbp_ideal = NULL;

  for (int i = 0; i < aintNumberOfThreadGroups; i++) 
  {
    if (ScoreBlkGroup[i]->kbp_ideal) 
    {
      lbsbpBlastScoreBlk->kbp_ideal = ScoreBlkGroup[i]->kbp_ideal;
    }
  }


  lbsbpBlastScoreBlk->ambig_size = ScoreBlkGroup[0]->ambig_size;
  lbsbpBlastScoreBlk->ambig_occupy = ScoreBlkGroup[0]->ambig_occupy;
  lbsbpBlastScoreBlk->ambiguous_res = ScoreBlkGroup[0]->ambiguous_res;
  lbsbpBlastScoreBlk->round_down = ScoreBlkGroup[0]->round_down;

  return lbsbpBlastScoreBlk;
}







BlastHSPResults *gbhrpfncCombineThreadGroupBlastHSPResults(int aintNumberOfThreadGroups, BlastHSPResults **abhrppBlastHSPResults) 
{
  int lintTotalQueries = 0;
             
  for (int i=0; i < aintNumberOfThreadGroups; i++) {
    //std::cout<<"HSPResults["<<i<<"]->num_queries= "<<abhrppBlastHSPResults[i]->num_queries<<std::endl;  std::cout.flush();
    lintTotalQueries += abhrppBlastHSPResults[i]->num_queries;
  }

  int lintCount = 0;
  //create new BlastHSPResults object and initialize with total number of queries we will be merging 
  BlastHSPResults *lbhrpBlastHSPResults = Blast_HSPResultsNew(lintTotalQueries);
  //std::cout<<"lbhrpBlastHSPResults->num_queries=  " <<lbhrpBlastHSPResults->num_queries <<std::endl;  std::cout.flush();
 
  int lintQueryOffset = 0;
  //loop over thread groups
  for (int ithgrp=0; ithgrp < aintNumberOfThreadGroups; ithgrp++) 
  {
    //loop over number of hitlist arrays
    int lintBlastHSPResultsQueryIndex;
    for (lintBlastHSPResultsQueryIndex = 0; 
	 lintBlastHSPResultsQueryIndex < (abhrppBlastHSPResults[ithgrp]->num_queries); 
	 ++lintBlastHSPResultsQueryIndex) 
    {
#if USE_SHALLOW_COPY
      //shallow copy hitlist_array
      //loop over hsp lists (arrays) and add offset for each thread group
      BlastHitList* old_hlst = abhrppBlastHSPResults[ithgrp]->hitlist_array[lintBlastHSPResultsQueryIndex];
      lbhrpBlastHSPResults->hitlist_array[lintBlastHSPResultsQueryIndex + lintCount] = old_hlst;
      BlastHitList* new_hlst = lbhrpBlastHSPResults->hitlist_array[lintBlastHSPResultsQueryIndex + lintCount];
      for (int ihsplst=0; ihsplst < new_hlst->hsplist_count; ihsplst++) 
      {
        new_hlst->hsplist_array[ihsplst]->query_index += lintQueryOffset;
      }
#else
      //deep copy hitlist array
      BlastHitList* old_hlst = abhrppBlastHSPResults[ithgrp]->hitlist_array[lintBlastHSPResultsQueryIndex];
      //std::cout<<"gbhrpfncCombineThreadGroupBlastHSPResults checkpt 3\n"<<std::cout.flush();
      //Allocate memory for new hitlist
      BlastHitList* new_hlst = Blast_HitListNew(old_hlst->hsplist_count);
      new_hlst->hsplist_count   = old_hlst->hsplist_count;
      new_hlst->hsplist_max     = old_hlst->hsplist_max;
      new_hlst->worst_evalue    = old_hlst->worst_evalue;
      new_hlst->low_score       = old_hlst->low_score;
      new_hlst->heapified       = old_hlst->heapified;
      new_hlst->hsplist_current = old_hlst->hsplist_current;
      //#if 1
      //works
      //allocate memory for hsplist array
      new_hlst->hsplist_array = (BlastHSPList**)calloc( old_hlst->hsplist_current, sizeof(BlastHSPList*) );
      for (int ihsplst=0; ihsplst < new_hlst->hsplist_count; ihsplst++) 
      {
	new_hlst->hsplist_array[ihsplst]= BlastHSPListDuplicate(old_hlst->hsplist_array[ihsplst]);
	//add offset to list query indices
        new_hlst->hsplist_array[ihsplst]->query_index += lintQueryOffset;
      }
      //std::cout<<"gbhrpfncCombineThreadGroupBlastHSPResults checkpt 8\n"<<std::cout.flush();
      //#else
      ////works
      //new_hlst->hsplist_array = old_hlst->hsplist_array;
      //for (int ihsplst=0; ihsplst < new_hlst->hsplist_count; ihsplst++) 
      //{
      //new_hlst->hsplist_array[ihsplst]->query_index += lintQueryOffset;
      //}
      //#endif
      //std::cout<<"gbhrpfncCombineThreadGroupBlastHSPResults checkpt 9\n"<<std::cout.flush();
      lbhrpBlastHSPResults->hitlist_array[lintBlastHSPResultsQueryIndex + lintCount] = new_hlst;
      //std::cout<<"gbhrpfncCombineThreadGroupBlastHSPResults checkpt 10\n"<<std::cout.flush();
#endif
    }
    
    //increment starting query count for next thread group
    lintCount += lintBlastHSPResultsQueryIndex;
    lintQueryOffset += lbhrpBlastHSPResults->num_queries;
  }

  return lbhrpBlastHSPResults; 
}


BlastQueryInfo *gbqipfncCombineThreadGroupBlastQueryInfo(int aintNumberOfThreadGroups, BlastQueryInfo **abqippBlastQueryInfo) 
{
  int lintTotalNumberOfQueries = 0;
  for (int i=0; i < aintNumberOfThreadGroups; i++) {
    lintTotalNumberOfQueries += abqippBlastQueryInfo[i]->num_queries;
  }

  //if restart and we are skipping this batch return an empty struct
  if(lintTotalNumberOfQueries==0){ return  BlastQueryInfoNew(eBlastTypeBlastx,0);}

  int lintCount = 0;

  BlastQueryInfo *lbqipBlastQueryInfo = BlastQueryInfoNew(eBlastTypeBlastx,lintTotalNumberOfQueries);

  lbqipBlastQueryInfo->first_context = 0;
  lbqipBlastQueryInfo->last_context = lintTotalNumberOfQueries-1;
  lbqipBlastQueryInfo->num_queries = lintTotalNumberOfQueries;
  lbqipBlastQueryInfo->max_length = abqippBlastQueryInfo[0]->max_length;
  lbqipBlastQueryInfo->pattern_info = abqippBlastQueryInfo[0]->pattern_info;

  int lintQueryOffset = 0;
  int lintQueryIndex = 0;

  //ceb the following code does a deep copy of each QueryInfo object
  for (int fintThreadGroupIndex = 0; fintThreadGroupIndex<aintNumberOfThreadGroups; fintThreadGroupIndex++) 
  {
    int fintBlastQueryInfoQueryIndex;
    //For each QueryInfo object loop over num_queries or contexts (same in this case)
    for (fintBlastQueryInfoQueryIndex = 0; 
	 fintBlastQueryInfoQueryIndex < (abqippBlastQueryInfo[fintThreadGroupIndex]->num_queries); 
	 ++fintBlastQueryInfoQueryIndex) 
    {
      //ceb record query length for all contexts
      lbqipBlastQueryInfo->contexts[fintBlastQueryInfoQueryIndex + lintCount].query_length =
	abqippBlastQueryInfo[fintThreadGroupIndex]->contexts[fintBlastQueryInfoQueryIndex].query_length;//ceb

      //record new offset for combined context list
      lbqipBlastQueryInfo->contexts[fintBlastQueryInfoQueryIndex + lintCount].query_offset = lintQueryOffset;

      //increment offset
      lintQueryOffset += abqippBlastQueryInfo[fintThreadGroupIndex]->contexts[fintBlastQueryInfoQueryIndex].query_length;

      //compute max length of all combined contexts
      if ( (abqippBlastQueryInfo[fintThreadGroupIndex]->contexts[fintBlastQueryInfoQueryIndex].query_length) > lbqipBlastQueryInfo->max_length ) 
      {
        lbqipBlastQueryInfo->max_length = 
	  abqippBlastQueryInfo[fintThreadGroupIndex]->contexts[fintBlastQueryInfoQueryIndex].query_length;
      }

      lbqipBlastQueryInfo->contexts[fintBlastQueryInfoQueryIndex + lintCount].eff_searchsp = 
	abqippBlastQueryInfo[fintThreadGroupIndex]->contexts[fintBlastQueryInfoQueryIndex].eff_searchsp;

      lbqipBlastQueryInfo->contexts[fintBlastQueryInfoQueryIndex + lintCount].length_adjustment = 
	abqippBlastQueryInfo[fintThreadGroupIndex]->contexts[fintBlastQueryInfoQueryIndex].length_adjustment;

      lbqipBlastQueryInfo->contexts[fintBlastQueryInfoQueryIndex + lintCount].query_index = lintQueryIndex;

      //increment combined context index
      lintQueryIndex++;
      //copy frame (int)
      lbqipBlastQueryInfo->contexts[fintBlastQueryInfoQueryIndex + lintCount].frame = 
	abqippBlastQueryInfo[fintThreadGroupIndex]->contexts[fintBlastQueryInfoQueryIndex].frame;

      lbqipBlastQueryInfo->contexts[fintBlastQueryInfoQueryIndex + lintCount].is_valid = 
	abqippBlastQueryInfo[fintThreadGroupIndex]->contexts[fintBlastQueryInfoQueryIndex].is_valid;
    }
    lintCount += fintBlastQueryInfoQueryIndex;
  }
  return lbqipBlastQueryInfo; 
}


CRef<CBlastQueryVector> gbqvfncCombineReplicationGroupQueries(int aintNumberOfThreadGroups,
							      char** buffer,
							      int gbl_bsize,
							      int* lintpBatchStartPosition,
							      int* lintpBatchEndPosition,
							      CBlastInputSourceConfig* iconfig,
							      int QueryBatchSize,
							      int gbl_query_count,
							      CRef<CScope> gbl_scope)
{
  //std::cout<<"gbqvfncCombineReplicationGroupQueries checkpt1 QueryBatchSize="<< QueryBatchSize<<"\n";std::cout.flush();
  //std::cout<<"gbqvfncCombineReplicationGroupQueries checkpt2 gbl_query_count="<<gbl_query_count<<"\n";std::cout.flush();

  //For tid==0 set up global query input stream
  
  //copy thread local query buffers into single global query buffer
  char*gbl_buffer = (char*) calloc( gbl_bsize+1, sizeof(char));
  int sumSize=0;
  //tack each lcl_buffer onto the gbl_buffer
  for ( int i=0; i < num_thread_groups; ++i ){
    int segSize = (lintpBatchEndPosition[i] - lintpBatchStartPosition[i]);
    //std::cout<<"lintpBatchEndPosition["<<i<<"]="<<lintpBatchEndPosition[i] <<" lintpBatchStartPosition["<<i<<"]"<<lintpBatchStartPosition[i]<<std::endl<<std::cout.flush();
    //std::cout<<"segSize= "<<segSize<<std::endl;std::cout.flush();
    sumSize+=segSize;
    strncat(gbl_buffer,&(buffer[i][ lintpBatchStartPosition[i]]),segSize);
    //only do the next step if we are not at the end of the query set
    if (num_thread_groups>1)
    {
      //tack on carriage return between local buffer segments
      strncat(gbl_buffer,"\n",1);sumSize+=1;
      //if(i==num_thread_groups-1){strncat(gbl_buffer,">",1);sumSize+=1;}
    }
    //char*tmp_buffer = (char*) calloc( gbl_bsize+1, sizeof(char));
    //strncat(tmp_buffer,&(buffer[i][ lintpBatchStartPosition[i]]),segSize);
    //std::cout<<"\n\nlcl_buffer["<<i<<"]\n"<<tmp_buffer<<std::endl<<std::endl;std::cout.flush();
  }

  //std::cout<<"sumSize =\n"<<sumSize<<std::endl;std::cout.flush();
  //std::cout<<"\n\ngbl_buffer\n"<< gbl_buffer << std::endl;std::cout.flush();

//ceb
//gbl_buffer has all expected queries
//by the time GetNextSeqBatch is called it is one short of exepected number
//This is causing the bug downstream

  //buffer[i] contains full query input from each file[i]
  membuf gbl_query_membuf(&(gbl_buffer[0]),&(gbl_buffer[sumSize-1]));

  auto_ptr<CNcbiIstream> gbl_query_input_stream;
  gbl_query_input_stream.reset(new CNcbiIstream(&gbl_query_membuf));
  
  if(IsIStreamEmpty(*gbl_query_input_stream)){ERR_POST(Warning << "Query is Empty!");}
  
  iconfig->SetLocalIdCounterInitValue(gbl_query_count+1);//ceb
  
  CBlastFastaInputSource gbl_fasta(*gbl_query_input_stream, *iconfig);
  
  CBlastInput gbl_input(&gbl_fasta, QueryBatchSize*num_thread_groups);
  //test to make sure we get all batches in gbl_input
  //gbl_input.SetBatchSize(QueryBatchSize*num_thread_groups*10);
  //CRef<CBlastQueryVector> gbl_query_batch(gbl_input.GetNextSeqBatch(*gbl_scope));
  CRef<CBlastQueryVector> gbl_query_batch(gbl_input.GetAllSeqs(*gbl_scope));

  //std::cout<<"gbqvfncCombineReplicationGroupQueries checkpt3 gbl_query_batch->Size()="<<gbl_query_batch->Size()<<"\n";std::cout.flush();
  //ceb test 
  //for(int i=0;i<gbl_query_batch->Size();++i)
  //{
  //  std::cout<<"gbl_query_batch["<<i<<"].GetQueryId() ="<<(((gbl_query_batch.GetObject())[i].GetObject()).GetQueryId().GetObject()).AsFastaString()<<"\n"<<std::cout.flush();
  //  //std::cout<<"gbl_query_batch["<<i<<"].GetLength() ="<<(((gbl_query_batch.GetObject())[i].GetObject()).GetLength())<<"\n"<<std::cout.flush();
															//}

  return gbl_query_batch;
};


//Function to create a search resluts set from the  4 combined objects

CRef<CSearchResultSet> gsrsfncCollectResults(BlastHSPResults* hsp_results,
					     BlastQueryInfo* m_QueryInfo,
					     BlastScoreBlk* m_ScoreBlk,
					     CRef<IQueryFactory> m_QueryFactory,
					     const CBlastOptions* m_Options,
					     CRef<IBlastSeqInfoSrc> m_SeqInfoSrc)
{
  //std::cout<<"gsrsfncCollectResults checkpt 1 \n";std::cout.flush();
    int hitlist_size_backup = m_Options->GetHitlistSize();
    EResultType m_ResultType = eDatabaseSearch;
    TSearchMessages m_Messages;// = ;//TSearchMessages&

    // This is the data resulting from the traceback phase (before it is converted to ASN.1).
    // We wrap it this way so it is released even if an exception is thrown below.

    // --shedsaw-- I think this may be problematic as it might call the free'ing function when
    // the object goes out of scope.

    //CRef< CStructWrapper<BlastHSPResults> > HspResults;
    //HspResults.Reset(WrapStruct(hsp_results, Blast_HSPResultsFree));

    _ASSERT(m_SeqInfoSrc);
    _ASSERT(m_QueryFactory);

//ceb causing segfault downstream
    CRef<ILocalQueryData> qdata = m_QueryFactory->MakeLocalQueryData(m_Options);



    m_SeqInfoSrc->GarbageCollect();
    vector<TSeqLocInfoVector> subj_masks;
    //converts hsp object into align vector
    //std::cout<<"gsrsfncCollectResults checkpt 5 \n";std::cout.flush();

    TSeqAlignVector aligns =
        LocalBlastResults2SeqAlign(hsp_results,
                                   *qdata,
                                   *m_SeqInfoSrc,
                                   m_Options->GetProgramType(),
                                   m_Options->GetGappedMode(),
                                   m_Options->GetOutOfFrameMode(),
                                   subj_masks,
                                   m_ResultType);
    //std::cout<<"gsrsfncCollectResults checkpt 6 \n";std::cout.flush();

    //query factory
    vector< CConstRef<CSeq_id> > query_ids;
    query_ids.reserve(aligns.size());
    //std::cout<<"gsrsfncCollectResults checkpt 7 \n";std::cout.flush();

    for (size_t i = 0; i < qdata->GetNumQueries(); i++){
      //std::cout<<"gbl Q_id: "<<qdata->GetSeq_loc(i)->GetId()<<std::endl;
        query_ids.push_back(CConstRef<CSeq_id>(qdata->GetSeq_loc(i)->GetId()));
    }
    //std::cout<<"gsrsfncCollectResults checkpt 8 \n";std::cout.flush();

    return BlastBuildSearchResultSet(query_ids,
                                     m_ScoreBlk,//->GetPointer(),
                                     m_QueryInfo,
                                     m_Options->GetProgramType(),
                                     aligns,
                                     m_Messages,
                                     subj_masks,
                                     NULL,
                                     m_ResultType);
}


#endif /* SKIP_DOXYGEN_PROCESSING */
