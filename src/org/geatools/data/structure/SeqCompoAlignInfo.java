package org.geatools.data.structure;
import java.util.Comparator;
/*******************************************************************************
 *  ========================================================================
 *  GEATools : a free Genomic Event Analysis Tools for the Java(tm) platform
 *  ========================================================================
 *
 *  (C) Copyright 2016, by Qi Wang and Contributors.
 *
 *  This file is part of GEATools.
 *
 *  GEATools is free software: you can redistribute it and/or modify it
 *  under the terms of the GNU Lesser General Public License as published by
 *  the Free Software Foundation; either version 2.1 of the License, or
 *  (at your option) any later version.
 *
 *  This library is distributed in the hope that it will be useful, but
 *  WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
 *  or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
 *  License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>
 *
 *  [Oracle and Java are registered trademarks of Oracle and/or its affiliates. 
 *  Other names may be trademarks of their respective owners.]
 *
 *  (C) Copyright 2000-2014, by Original Author and Contributors.
 *
 *  Original Author: Qi Wang;
 *  Contributor(s): 
 *  Changes (from 1-Jan-2016)
 *  ---------------------------------------------------------------------------
 *******************************************************************************/
import java.util.List;

public class SeqCompoAlignInfo extends SeqInfo{
	  public int seqIndex;  
	  public List<Integer> seqAlignSStart;
	  public List<Integer> seqAlignSEnd;
	  public List<Integer> seqAlignQStart;
	  public List<Integer> seqAlignQEnd;
	  public List<Boolean> isExactMatch;
	  public int baitBrkSPos=-1;
	  public int baitBrkQPos=-1;
	  /*
	  public int genomeStart=-1;
	  public int genomeEnd=-1;
	  public float genomeAlignScore;
	  public String chr="";
	  public String chrStrand=""; 
	  */ 
	  
	  public static class CompSeqEncode implements Comparator<SeqCompoAlignInfo> {
			
			private int mod = 1;
			public CompSeqEncode(boolean desc) {
			  if (desc) mod =-1;
			}
			        
			public int compare(SeqCompoAlignInfo  seq1, SeqCompoAlignInfo  seq2){
			  return  mod*Long.valueOf(seq1.seqNumEncode).compareTo(seq2.seqNumEncode);
			}
			
	  }
}
