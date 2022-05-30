package org.geatools.data.structure;
import java.util.ArrayList;
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

public class SeqRocket implements Cloneable{
	public String rocketName="";
	public String expName="";
	public boolean isActive=true;
	public boolean isDone=false;	
	public List<SeqCompoRecognizer> seqRecognizers;
	public SeqCompoFeatures seqCompoFeatures;
	public List<String> savedCompoAlignedSeqFiles;
	public String finalSeqCompo;
	public String outDir;
	public String recognizedSeqFileMasked;
	public String recognizedSeqFileTrimmed;
	public String recognizedSeqFile;
	public String recognizedSeqHTMLFile;
	public String minLeftRecognizerSeq;
	public String allRightRecognizerSeq;
	public int recognizedSeqCount=0;
	public int exactSeqCount=0;
	public int nonExactSeqCount=0;
	public double recognizedSeqRPM=0;
	public boolean saveRecognizedSeq=true;
	public boolean saveRecognizedSeqAsFASTA=true;
	public boolean saveRecognizedSeqAsFASTQ=false;
	public boolean saveRecognizedSeqAsHTML;
	public boolean isRecognized=false;
	public boolean isEndTrimmed=false;
	public boolean isEndMasked=false;
	
	public void setSeqRecognizer(List<SeqCompoRecognizer> seqRecognizerList){
		this.seqRecognizers=seqRecognizerList;
	}
	public void setSeqCompoType(SeqCompoFeatures seqCompoFeatureList){
		this.seqCompoFeatures=seqCompoFeatureList;
	}	

	
	@Override
	public Object clone() throws CloneNotSupportedException {
		SeqRocket cloned = (SeqRocket)super.clone();
		// deep copy of object field
		List<SeqCompoRecognizer> clonedSeqRecognizer=new ArrayList<SeqCompoRecognizer>();
		for(SeqCompoRecognizer seqRecognizer:cloned.seqRecognizers){
		  SeqCompoRecognizer recognizer=(SeqCompoRecognizer) seqRecognizer.clone();
		  clonedSeqRecognizer.add(recognizer);
		}
		cloned.setSeqRecognizer(clonedSeqRecognizer);
		
		SeqCompoFeatures seqCompo=(SeqCompoFeatures) cloned.seqCompoFeatures.clone();
		cloned.setSeqCompoType(seqCompo);
		
	    return cloned;
	}
}


