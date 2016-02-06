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
package org.geatools.seqprocess;

import java.util.ArrayList;
import java.util.List;

import org.geatools.data.structure.SeqInfo;
import org.geatools.operation.FileOperate;


public class SeqDupFilter{

 public SeqDupFilter(){

 }
 
 public static int filterExactDup(String seqFile,String subSeqNameFile,int seqNameColIdx){
	
     int num=0;
	 List<SeqInfo> subSeqObjList
	     =SeqOperation.extractSubSeq(seqFile,subSeqNameFile,seqNameColIdx);
	
	 List<SeqInfo> seqObjList = SeqOperation.getSeqObj(seqFile);	
	 List<String> excludedSeqNameList=new ArrayList<String>();
	 List<SeqInfo> dupSeqObjList=new ArrayList<SeqInfo>(); 
	 for(SeqInfo subSeqObj:subSeqObjList){
		for(SeqInfo seqObj:seqObjList){				
			if(seqObj.seq.equalsIgnoreCase(subSeqObj.seq)){
			  dupSeqObjList.add(seqObj);
			  excludedSeqNameList.add(seqObj.seqName);
			}
		}
	 }
	 subSeqObjList=null;
	 String format=FileOperate.getFileFormat(seqFile);
	 String outSeqFile=subSeqNameFile+"Dup."+format;
	 SeqOperation.saveSeqList(dupSeqObjList, outSeqFile);
	 dupSeqObjList=null;
	 
	 List<SeqInfo> filterSeqObjList=SeqOperation.excludeSeq(seqObjList,excludedSeqNameList);
	 outSeqFile=seqFile.substring(0,seqFile.lastIndexOf("."))+".clean."+format;
	 SeqOperation.saveSeqList(filterSeqObjList, outSeqFile);
	 
	 num=filterSeqObjList.size();
	 
	 filterSeqObjList=null;
	 
	 return num;
  }
 
 /*
  public static List<SeqInfo> filterExactDup(List<SeqInfo> seqObjList,
		  String subSeqNameFile,int seqNameColIdx){
		
	   
	 List<SeqInfo> subSeqObjList
	     =SeqOperation.extractSubSeq(seqObjList,subSeqNameFile,seqNameColIdx);
	
		
	 List<String> excludedSeqNameList=new ArrayList<String>();
	 List<SeqInfo> dupSeqObjList=new ArrayList<SeqInfo>(); 
	 for(SeqInfo subSeqObj:subSeqObjList){
		for(SeqInfo seqObj:seqObjList){				
			if(seqObj.seq.equalsIgnoreCase(subSeqObj.seq)){
				dupSeqObjList.add(seqObj);
				excludedSeqNameList.add(seqObj.seqName);
			}
		}
	 }
	 subSeqObjList=null;
	 dupSeqObjList=null;
	 
	 List<SeqInfo> filterSeqObjList=SeqOperation.excludeSeq(seqObjList,excludedSeqNameList);
	 //outSeqFile=seqFile.substring(0,seqFile.lastIndexOf("."))+".clean."+format;
	 //SeqOperation.saveSeqList(subSeqObjList, outSeqFile);
	 
	 return filterSeqObjList;
  }
  */
}