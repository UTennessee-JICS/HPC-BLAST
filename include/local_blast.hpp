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


/* $Id: local_blast.hpp 396870 2013-04-23 14:03:03Z madden $
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
 * Author: Christiam Camacho, Kevin Bealer
 *
 */

/** @file local_blast.hpp
 * Main class to perform a BLAST search on the local machine.
 */

#ifndef ALGO_BLAST_API___LOCAL_BLAST_HPP
#define ALGO_BLAST_API___LOCAL_BLAST_HPP

#include <algo/blast/api/prelim_stage.hpp>
#include <algo/blast/api/traceback_stage.hpp>
#include <algo/blast/api/blast_seqinfosrc.hpp>

/** @addtogroup AlgoBlast
 *
 * @{
 */
class CBlastFilterTest;

BEGIN_NCBI_SCOPE
BEGIN_SCOPE(blast)

/// Returns the optimal chunk size for a given task
/// @param program BLAST task [in]
NCBI_XBLAST_EXPORT
size_t
SplitQuery_GetChunkSize(EProgram program);

/// Class to perform a BLAST search on local BLAST databases
/// Note that PHI-BLAST can be run using this class also, one only need to
/// configure it as a regular blastp or blastn search and set the pattern in
/// the CBlastOptionsHandle object
/// @todo should RPS-BLAST be moved out of this object?
class NCBI_XBLAST_EXPORT CLocalBlast : public CObject, public CThreadable
{
public:
    /// Constructor with database description
    /// @param query_factory query sequence(s) [in]
    /// @param opts_handle BLAST options handle [in]
    /// @param dbinfo description of BLAST database to search [in]
    CLocalBlast(CRef<IQueryFactory> query_factory,
                CRef<CBlastOptionsHandle> opts_handle,
                const CSearchDatabase& dbinfo);

    /// Constructor with subject adapter (@sa CLocalDbAdapter)
    /// @param query_factory query sequence(s) [in]
    /// @param opts_handle BLAST options handle [in]
    /// @param db subject adapter [in]
    CLocalBlast(CRef<IQueryFactory> query_factory,
                CRef<CBlastOptionsHandle> opts_handle,
                CRef<CLocalDbAdapter> db);

    /// Constructor with database description
    /// @param query_factory query sequence(s) [in]
    /// @param opts_handle BLAST options handle [in]
    /// @param seqsrc BlastSeqSrc object to search [in]
    /// @param seqInfoSrc user-specified IBlastSeqInfoSrc [in]
    CLocalBlast(CRef<IQueryFactory> query_factory,
                CRef<CBlastOptionsHandle> opts_handle,
                BlastSeqSrc* seqsrc,
                CRef<IBlastSeqInfoSrc> seqInfoSrc);
    
    /// Executes the search
    CRef<CSearchResultSet> Run();
    
    /// Set a function callback to be invoked by the CORE of BLAST to allow
    /// interrupting a BLAST search in progress.
    /// @param fnptr pointer to callback function [in]
    /// @param user_data user data to be attached to SBlastProgress structure 
    /// [in]
    /// @return the previously set TInterruptFnPtr (NULL if none was 
    /// provided before)
    TInterruptFnPtr SetInterruptCallback(TInterruptFnPtr fnptr,
                                         void* user_data = NULL) {
        _ASSERT(m_PrelimSearch);
        return m_PrelimSearch->SetInterruptCallback(fnptr, user_data);
    }
  
    /// Retrieve any error/warning messages that occurred during the search
    TSearchMessages GetSearchMessages() const;

    /// Retrieve the number of extensions performed during the search
    Int4 GetNumExtensions();

    /// Get the diagnostics structure (deep copy, needs to be deleted by caller)
    BlastDiagnostics* GetDiagnostics();

    //ceb
    char* GetHSPBuffer();
    char* GetQueryInfoBuffer();
    char* GetScoreBlkBuffer();

    size_t GetHSPBSize();
    size_t GetQueryInfoBSize();
    size_t GetScoreBlkBSize();

    //ceb New object to save for thread level parallelism
    //CRef< CStructWrapper<BlastHSPResults> > GetHSPResults();
    BlastHSPResults* GetHSPResults();
    BlastQueryInfo*  GetQueryInfo();
    BlastScoreBlk*   GetScoreBlk();


private:
    /// Query factory from which to obtain the query sequence data
    CRef<IQueryFactory> m_QueryFactory;
    
    /// Options to use
    CRef<CBlastOptions> m_Opts;

    /// Internal core data structures which are used in the preliminary and
    /// traceback stages of the search
    CRef<SInternalData> m_InternalData;

    /// Object which runs the preliminary stage of the search
    CRef<CBlastPrelimSearch> m_PrelimSearch;

    /// Object which runs the traceback stage of the search
    CRef<CBlastTracebackSearch> m_TbackSearch;
    
    /// Local DB adaptor (if one was) passed to constructor.
    CRef<CLocalDbAdapter> m_LocalDbAdapter;

    /// User-specified IBlastSeqInfoSrc implementation
    /// (may be used for non-standard databases, etc.)
    CRef<IBlastSeqInfoSrc> m_SeqInfoSrc;

    /// Warnings and error messages
    TSearchMessages                 m_Messages;

    //ceb Buffers for the data structures.
    //we may be able to remove these buffers by moving the
    // packing and unpacking into the hpc_blastp_app
    char * tHSPResultsBuffer;
    char * tQueryInfoBuffer;
    char * tBlastScoreBlkBuffer;
    //ceb Size of the buffers.
    size_t tHSPResultsBSize;
    size_t tQueryInfoBSize;
    size_t tBlastScoreBlkBSize;

    //ceb
    //Handles to pass out objects for thread level parallelism
    BlastHSPResults* HSPResults;
    //CRef< CStructWrapper<BlastHSPResults> > HSPResults;
    BlastQueryInfo* QueryInfo;
    BlastScoreBlk* ScoreBlk;

    friend class ::CBlastFilterTest;
    friend class CBl2Seq;
};

inline TSearchMessages
CLocalBlast::GetSearchMessages() const
{
    return m_Messages;
}

inline char* 
CLocalBlast::GetHSPBuffer()
{
  return tHSPResultsBuffer;
}

inline size_t 
CLocalBlast::GetHSPBSize()
{
  return tHSPResultsBSize;
}

inline char* 
CLocalBlast::GetQueryInfoBuffer()
{
  return tQueryInfoBuffer;
}

inline size_t 
CLocalBlast::GetQueryInfoBSize()
{
  return tQueryInfoBSize;
}

inline char* 
CLocalBlast::GetScoreBlkBuffer()
{
  return tBlastScoreBlkBuffer;
}

inline size_t 
CLocalBlast::GetScoreBlkBSize()
{
  return tBlastScoreBlkBSize;
}

//ceb New object to save for thread level parallelism
//inline CRef< CStructWrapper<BlastHSPResults> > 
inline BlastHSPResults* 
CLocalBlast::GetHSPResults()
{
   return HSPResults;
}
inline BlastQueryInfo* 
CLocalBlast::GetQueryInfo()
{
   return QueryInfo;
}
inline BlastScoreBlk* 
CLocalBlast::GetScoreBlk()
{
   return ScoreBlk;
}

END_SCOPE(BLAST)
END_NCBI_SCOPE

/* @} */

#endif /* ALGO_BLAST_API___LOCAL_BLAST__HPP */
