package org.geatools.data.structure;

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
import java.util.Comparator;

public class SeqInfo{
	  public int seqID;
	  public String seqName;
	  public String seqIdentifier;
	  public int seqLength;
	  public String seq; 
	  public String seqQualityEncode; 
	  public long seqNumEncode=0;
	  public long seqNumRevEncode=0;
	  public boolean isDup=false;
	  public int dupNum=-1;
	  
	  public static class CompSeqEncode implements Comparator<SeqInfo> {
			
			private int mod = 1;
			public CompSeqEncode(boolean desc) {
			  if (desc) mod =-1;
			}
			        
			public int compare(SeqInfo  seq1, SeqInfo  seq2){
			  return  mod*Long.valueOf(seq1.seqNumEncode).compareTo(seq2.seqNumEncode);
			}
			
	  }
	  
	  public static class CompSeqRevEncode implements Comparator<SeqInfo> {
			
			private int mod = 1;
			public CompSeqRevEncode(boolean desc) {
			  if (desc) mod =-1;
			}
			        
			public int compare(SeqInfo  seq1, SeqInfo  seq2){
			  return  mod*Long.valueOf(seq1.seqNumRevEncode).compareTo(seq2.seqNumRevEncode);
			}
			
	  }

	  
	  public static class CompSeqID implements Comparator<SeqInfo> {
			
			private int mod = 1;
			public CompSeqID(boolean desc) {
			  if (desc) mod =-1;
			}
			        
			public int compare(SeqInfo  seq1, SeqInfo  seq2){
			  return  mod*Integer.valueOf(seq1.seqID).compareTo(seq2.seqID);
			}
			
	  }

}