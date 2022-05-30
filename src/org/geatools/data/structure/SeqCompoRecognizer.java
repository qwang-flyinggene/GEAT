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
public class SeqCompoRecognizer extends SeqInfo implements Cloneable{
	  public int index=0;
	  public String rawSeq;
	  public int rawSeqLength;
	  public String side=null;
	  public String seqType;
	  public String seqFASTAFile;	 
	  public String exactAlignedSeqFile;	  
	  public String alignedSeqFile;
	  public String alignedSeqFileTrimmed;
	  public String alignedSeqFileMasked;
	  public String alignedSeqHTMLFile;
	  public int exactMaxStart=1;
	  public int exactMinStart=1;
	  public int minSeqLen=8;
	  public int territoryLen=0;
	  public float territoryPercent=0.6f;
	  public double minAlignRatio=0.0d;
	  public double maxMismatchRatio=0.0d;
	  public double maxGapRatio=0.0d;
	  public int minAlignLen=0;
	  public int maxMismatchNum=0;
	  public int maxGapNum=0;  
	  public int maxQGapExtension=0;
	  public int maxSGapExtension=0;
	  public int maxQStart=0;	  
	  public int minQStart=0;
	  public int maxQEnd=0;	  
	  public int minQEnd=0;
	  public int maxSStart=0;
	  public int minSStart=0;
	  public int maxSEnd=0;
	  public int minSEnd=0;
	  public int blastWordSize=7; 
	  public String blastTask="blastn-short";  
	  public boolean leftSubForBLAST=false;
	  public boolean leftShiftForNext=false;
	  public int trimLeftShift=0;
	  public boolean rightSubForBLAST=false;
	  public boolean rightShiftForNext=false;	
	  public int trimRightShift=0;	  
	  public boolean saveAsFASTA=false;
	  public boolean saveAsFASTQ=false;
	  public boolean isFASTASaved=false;
	  public boolean isFASTQSaved=false;
	  public boolean leftTrimSave=false;
	  public boolean rightTrimSave=false;
	  public boolean leftMaskSave=false;
	  public boolean rightMaskSave=false;
	  public boolean done=false;
	  
	  @Override
	  public Object clone() throws CloneNotSupportedException {
			SeqCompoRecognizer cloned = (SeqCompoRecognizer)super.clone();
		   
			return cloned;
	  }
	    
}
