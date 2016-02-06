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

public class SeqChrAlignSite extends ChrSite{
	 
	 public int qSeqID;
	 public String qName;
	 public String sName;
	 public float identity;
	 public int alignLen;
	 public int mismatchNum;
	 public int gapNum;
	 public int qStart;
	 public int qEnd;
	 public int sStart;
	 public int sEnd;
	 public double eValue;
	 public double score;
	    
	 public static class CompScore implements Comparator<SeqChrAlignSite> {
		
		private int mod = 1;
		public CompScore(boolean desc) {
		  if (desc) mod =-1;
		}
		        
		public int compare(SeqChrAlignSite  site1, SeqChrAlignSite  site2){
		  return  mod*Double.valueOf(site1.score).compareTo(site2.score);
		}
		
	 }
	 
	 public static class CompMismatch implements Comparator<SeqChrAlignSite> {
			
		private int mod = 1;
		public CompMismatch(boolean desc) {
		  if (desc) mod =-1;
		}
		        
		public int compare(SeqChrAlignSite  site1, SeqChrAlignSite  site2){
		  return  mod*Integer.valueOf(site1.mismatchNum).compareTo(site2.mismatchNum);
		}
		
	 }

}
