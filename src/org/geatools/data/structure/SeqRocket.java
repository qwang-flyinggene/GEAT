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
	public boolean isActive=true;
	
	public SeqCompo seqComponent;
	public List<SeqCompoRecognizer> seqRecognizer;
	public SeqCompoType seqTypeInfo;
	public String barcodeRecognizedSeqFile;
	public String recognizedMaskSeqFile;
	public String recognizedTrimSeqFile;
	public String recognizedSeqFile;
	public int recognizedSeq;
	public int exactAlignedSeq=0;
	
	public void setSeqRecognizer(List<SeqCompoRecognizer> seqRecognizerList){
		this.seqRecognizer=seqRecognizerList;
	}
	
	public void setSeqComponent(SeqCompo seqCompo){
		this.seqComponent=seqCompo;
	}
	
	public void setSeqTypeInfo(SeqCompoType seqType){
		this.seqTypeInfo=seqType;
	}
	
	
	@Override
	public Object clone() throws CloneNotSupportedException {
		SeqRocket cloned = (SeqRocket)super.clone();
		// deep copy of object field
		List<SeqCompoRecognizer> clonedSeqRecognizer=new ArrayList<SeqCompoRecognizer>();
		for(SeqCompoRecognizer seqRecognizer:cloned.seqRecognizer){
		  SeqCompoRecognizer recognizer=(SeqCompoRecognizer) seqRecognizer.clone();
		  clonedSeqRecognizer.add(recognizer);
		}
		cloned.setSeqRecognizer(clonedSeqRecognizer);
		
		SeqCompo seqCompo=(SeqCompo) cloned.seqComponent.clone();
		cloned.setSeqComponent(seqCompo);
		
		SeqCompoType seqType=(SeqCompoType) cloned.seqTypeInfo.clone();
		cloned.setSeqTypeInfo(seqType);
		
	    return cloned;
	}
}


