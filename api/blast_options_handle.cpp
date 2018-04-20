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


/*  $Id: blast_options_handle.cpp 442346 2014-08-01 16:28:16Z mcelhany $
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

/// @file blast_options_handle.cpp
/// Implementation for the CBlastOptionsHandle and the
/// CBlastOptionsFactory classes.

#include <ncbi_pch.hpp>
#include <algo/blast/api/blast_exception.hpp>
#include <algo/blast/api/blast_options_handle.hpp>
#include <algo/blast/api/blast_prot_options.hpp>
#include <algo/blast/api/blastx_options.hpp>
#include <algo/blast/api/tblastn_options.hpp>
#include <algo/blast/api/rpstblastn_options.hpp>
#include <algo/blast/api/tblastx_options.hpp>
#include <algo/blast/api/blast_nucl_options.hpp>
#include <algo/blast/api/disc_nucl_options.hpp>
#include <algo/blast/api/blast_rps_options.hpp>
#include <algo/blast/api/blast_advprot_options.hpp>
#include <algo/blast/api/psiblast_options.hpp>
#include <algo/blast/api/phiblast_nucl_options.hpp>
#include <algo/blast/api/phiblast_prot_options.hpp>
#include <algo/blast/api/deltablast_options.hpp>

/** @addtogroup AlgoBlast
 *
 * @{
 */

BEGIN_NCBI_SCOPE
BEGIN_SCOPE(blast)

CBlastOptionsHandle::CBlastOptionsHandle(EAPILocality locality)
    : m_DefaultsMode(false)
{
    m_Opts.Reset(new CBlastOptions(locality));
}

CBlastOptionsHandle::CBlastOptionsHandle(CRef<CBlastOptions> opt)
    : m_Opts(opt),
      m_DefaultsMode(false)
{
}

void
CBlastOptionsHandle::SetDefaults()
{
    if (m_Opts->GetLocality() != CBlastOptions::eRemote) {
        m_Opts->SetDefaultsMode(true);
        SetLookupTableDefaults();
        SetQueryOptionDefaults();
        SetInitialWordOptionsDefaults();
        SetGappedExtensionDefaults();
        SetScoringOptionsDefaults();
        SetHitSavingOptionsDefaults();
        SetEffectiveLengthsOptionsDefaults();
        SetSubjectSequenceOptionsDefaults();
        m_Opts->SetDefaultsMode(false);
    }
    SetRemoteProgramAndService_Blast3();
}

bool
CBlastOptionsHandle::Validate() const
{
    return m_Opts->Validate();
}

char*
CBlastOptionsHandle::GetFilterString() const
{ 
    return m_Opts->GetFilterString(); /* NCBI_FAKE_WARNING */
}

void 
CBlastOptionsHandle::SetFilterString(const char* f, bool clear /* = true */) 
{
    m_Opts->SetFilterString(f, clear); /* NCBI_FAKE_WARNING */
}

CBlastOptionsHandle*
CBlastOptionsFactory::Create(EProgram program, EAPILocality locality)
{
    CBlastOptionsHandle* retval = NULL;

    switch (program) {
    case eBlastn:
	{
		CBlastNucleotideOptionsHandle* opts = 
			new CBlastNucleotideOptionsHandle(locality);
		opts->SetTraditionalBlastnDefaults();
		retval = opts;
        break;
	}

    case eBlastp:
        retval = new CBlastAdvancedProteinOptionsHandle(locality);
        break;

    case eBlastx:
        retval = new CBlastxOptionsHandle(locality);
        break;

    case eTblastn:
        retval = new CTBlastnOptionsHandle(locality);
        break;

    case eTblastx:
        retval = new CTBlastxOptionsHandle(locality);
        break;

    case eMegablast:
	{
		CBlastNucleotideOptionsHandle* opts = 
			new CBlastNucleotideOptionsHandle(locality);
		opts->SetTraditionalMegablastDefaults();
		retval = opts;
        break;
	}

    case eDiscMegablast:
        retval = new CDiscNucleotideOptionsHandle(locality);
        break;

    case eRPSBlast:
        retval = new CBlastRPSOptionsHandle(locality);
        break;

    case eRPSTblastn:
        retval = new CRPSTBlastnOptionsHandle(locality);
        break;
        
    case ePSIBlast:
        retval = new CPSIBlastOptionsHandle(locality);
        break;

    case ePSITblastn:
        retval = new CPSIBlastOptionsHandle(locality);
        (dynamic_cast<CPSIBlastOptionsHandle *> (retval))->SetPSITblastnDefaults();
        break;

    case ePHIBlastp:
        retval = new CPHIBlastProtOptionsHandle(locality);
        break;        

    case ePHIBlastn:
        retval = new CPHIBlastNuclOptionsHandle(locality);
        break;

    case eDeltaBlast:
        retval = new CDeltaBlastOptionsHandle(locality);
        break;

    case eVecScreen:
	{
		CBlastNucleotideOptionsHandle* opts = 
			new CBlastNucleotideOptionsHandle(locality);
		opts->SetVecScreenDefaults();
		retval = opts;
        break;
	}

    case eBlastNotSet:
        NCBI_THROW(CBlastException, eInvalidArgument,
                   "eBlastNotSet may not be used as argument");
        break;

    default:
        abort();    // should never happen
    }

    return retval;
}

set<string>
CBlastOptionsFactory::GetTasks(ETaskSets choice /* = eAll */)
{
    set<string> retval;
    if (choice == eNuclNucl || choice == eAll) {
        retval.insert("blastn");
        retval.insert("blastn-short");
        retval.insert("megablast");
        retval.insert("dc-megablast");
        retval.insert("vecscreen");
        // -RMH-
        retval.insert("rmblastn");
    }

    if (choice == eProtProt || choice == eAll) {
        retval.insert("blastp");
        retval.insert("blastp-short");
        retval.insert("blastp-fast");
    }

    if (choice == eAll) {
        retval.insert("psiblast");
        //retval.insert("phiblastn"); // not supported yet
        retval.insert("phiblastp");
        retval.insert("rpsblast");
        retval.insert("rpstblastn");
        retval.insert("blastx");
        retval.insert("blastx-fast");
        retval.insert("deltablast");
        retval.insert("tblastn");
        retval.insert("tblastn-fast");
        retval.insert("psitblastn");
        retval.insert("tblastx");
    }

    return retval;
}

string
CBlastOptionsFactory::GetDocumentation(const string& task_name)
{
    string task(task_name);
    NStr::ToLower(task);
    string retval;

    if (task == "blastn") {
        retval.assign("Traditional BLASTN requiring an exact match of 11");
    } else if (task == "blastn-short") {
        retval.assign("BLASTN program optimized for sequences shorter than ");
        retval += "50 bases";
    } else if (task == "vecscreen") {
        retval.assign("BLASTN with several options re-set for running VecScreen");
    } else if (task == "rmblastn") {
        retval.assign("BLASTN with complexity adjusted scoring and masklevel");
        retval += "filtering";
    } else if (task == "blastp") {
        retval.assign("Traditional BLASTP to compare a protein query to a ");
        retval += "protein database";
    } else if (task == "blastp-short") {
        retval.assign("BLASTP optimized for queries shorter than 30 residues");
    } else if (task == "blastp-fast") {
        retval.assign("BLASTP optimized for faster runtime");
    } else if (task == "blastx") {
        retval.assign("Search of a (translated) nucleotide query against a ");
        retval += "protein database";
    } else if (task == "blastx-fast") {
        retval.assign("Search of a (translated) nucleotide query against a ");
        retval += "protein database with parameters optimized for faster runtime";
    } else if (task == "dc-megablast") {
        retval.assign("Discontiguous megablast used to find more distant ");
        retval += "(e.g., interspecies) sequences";
    } else if (task == "megablast") {
        retval.assign("Traditional megablast used to find very similar ");
        retval += "(e.g., intraspecies or closely related species) sequences";
    } else if (NStr::StartsWith(task, "phiblast")) {
        retval.assign("Limits BLASTP search to those subjects with a ");
        retval += "pattern matching one in the query";
    } else if (task == "psiblast") {
        retval.assign("PSIBLAST that searches a (protein) profile against ");
        retval += "a protein database";
    } else if (task == "rpsblast") {
        retval.assign("Search of a protein query against a database of motifs");
    } else if (task == "rpstblastn") {
        retval.assign("Search of a (translated) nucleotide query against ");
        retval.append("a database of motifs");
    } else if (task == "tblastn") {
        retval.assign("Search of a protein query against a (translated) ");
        retval += "nucleotide database";
    } else if (task == "tblastn-fast") {
        retval.assign("Search of a protein query against a (translated) ");
        retval += "nucleotide database with parameters optimized for faster runtime";
    } else if (task == "psitblastn") {
        retval.assign("Search of a PSSM against a (translated) ");
        retval += "nucleotide database";
    } else if (task == "tblastx") {
        retval.assign("Search of a (translated) nucleotide query against ");
        retval += "a (translated) nucleotide database";
    } else if (task == "deltablast") {
        retval.assign("DELTA-BLAST builds profile using conserved domain ");
        retval += "and uses this profile to search protein database";
    } else {
        retval.assign("Unknown task");
    }
    return retval;
}

CBlastOptionsHandle*
CBlastOptionsFactory::CreateTask(string task, EAPILocality locality)
{
    CBlastOptionsHandle* retval = NULL;

    string lc_task(NStr::ToLower(task));
    //ThrowIfInvalidTask(lc_task);

    if (!NStr::CompareNocase(task, "blastn") || 
        !NStr::CompareNocase(task, "blastn-short") ||
        // -RMH-
        !NStr::CompareNocase(task, "rmblastn") ||
        !NStr::CompareNocase(task, "vecscreen"))
    {
        CBlastNucleotideOptionsHandle* opts = 
             dynamic_cast<CBlastNucleotideOptionsHandle*>
                (CBlastOptionsFactory::Create(eBlastn, locality));
        _ASSERT(opts);
        if (!NStr::CompareNocase(task, "blastn-short"))
        {
             opts->SetMatchReward(1);
             opts->SetMismatchPenalty(-3);
             opts->SetEvalueThreshold(50);
             opts->SetWordSize(7);
             opts->ClearFilterOptions();
        }
        else if (!NStr::CompareNocase(task, "vecscreen"))
        {
            opts->SetVecScreenDefaults();
        }else if ( !NStr::CompareNocase(task, "rmblastn") )
        {
            // -RMH- This blastn only supports full matrix scoring.
            opts->SetMatchReward(0);
            opts->SetMismatchPenalty(0);
        }
        retval = opts;
    }
    else if (!NStr::CompareNocase(task, "megablast"))
    {
         retval = CBlastOptionsFactory::Create(eMegablast, locality);
    }
    else if (!NStr::CompareNocase(task, "dc-megablast"))
    {
         retval = CBlastOptionsFactory::Create(eDiscMegablast, locality);
    }
    else if (!NStr::CompareNocase(task, "blastp") || 
             !NStr::CompareNocase(task, "blastp-short") ||
             !NStr::CompareNocase(task, "blastp-fast"))
    {
         CBlastAdvancedProteinOptionsHandle* opts =
               dynamic_cast<CBlastAdvancedProteinOptionsHandle*> 
                (CBlastOptionsFactory::Create(eBlastp, locality));
         if (task == "blastp-short") {
            opts->SetMatrixName("PAM30");
            opts->SetGapOpeningCost(9);
            opts->SetGapExtensionCost(1);
            opts->SetEvalueThreshold(20000);
            opts->SetWordSize(2);
            opts->ClearFilterOptions();
         } else if (task == "blastp-fast") {
            opts->SetWordSize(6);
            opts->SetOptions().SetLookupTableType(eCompressedAaLookupTable);
            opts->SetWordThreshold(21.0);
         }
         retval = opts;
    }
    else if (!NStr::CompareNocase(task, "psiblast"))
    {
         retval = CBlastOptionsFactory::Create(ePSIBlast, locality);
    }
    else if (!NStr::CompareNocase(task, "psitblastn"))
    {
         retval = CBlastOptionsFactory::Create(ePSITblastn, locality);
    }
    else if (!NStr::CompareNocase(task, "phiblastp"))
    {
         retval = CBlastOptionsFactory::Create(ePHIBlastp, locality);
    }
    else if (!NStr::CompareNocase(task, "rpsblast"))
    {
         retval = CBlastOptionsFactory::Create(eRPSBlast, locality);
    }
    else if (!NStr::CompareNocase(task, "rpstblastn"))
    {
         retval = CBlastOptionsFactory::Create(eRPSTblastn, locality);
    }
    else if (!NStr::CompareNocase(task, "blastx") ||
             !NStr::CompareNocase(task, "blastx-fast"))
    {
        CBlastxOptionsHandle* opts =
            dynamic_cast<CBlastxOptionsHandle*>
            (CBlastOptionsFactory::Create(eBlastx, locality));
        if (task == "blastx-fast") {
            opts->SetWordSize(6);
            opts->SetOptions().SetLookupTableType(eCompressedAaLookupTable);
            opts->SetWordThreshold(21.0);
        }
        retval = opts;
    }
    else if (!NStr::CompareNocase(task, "tblastn") ||
    		 !NStr::CompareNocase(task, "tblastn-fast"))
    {
    	CTBlastnOptionsHandle* opts =
    	            dynamic_cast<CTBlastnOptionsHandle*>
         (CBlastOptionsFactory::Create(eTblastn, locality));
    	if(task == "tblastn-fast") {
    		opts->SetWordSize(6);
            opts->SetOptions().SetLookupTableType(eCompressedAaLookupTable);
            opts->SetWordThreshold(21.0);
    	}
    	retval = opts;
    }
    else if (!NStr::CompareNocase(task, "tblastx"))
    {
         retval = CBlastOptionsFactory::Create(eTblastx, locality);
    }
    else if (!NStr::CompareNocase(task, "deltablast"))
    {
         retval = CBlastOptionsFactory::Create(eDeltaBlast, locality);
    }
    else
    {
        abort();    // should never get here
    }
    return retval;
}

END_SCOPE(blast)
END_NCBI_SCOPE


/* @} */
