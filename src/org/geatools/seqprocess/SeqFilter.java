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

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Calendar;
import java.util.Collections;
import java.util.List;

import org.geatools.GEAT;
import org.geatools.data.structure.SeqEncode;
import org.geatools.data.structure.SeqInfo;
import org.geatools.operation.FileOperation;
import org.geatools.operation.SeqOperation;


public class SeqFilter{
 
   static int maxEncodeBaseLen=19;
   static int seqSplitStep=1000000;
   static int encodeIntervalNum=20;
 
   public SeqFilter(){

   }
   
   public void setMaxEncodeBaseLen(int len){ 
	  if(len>19) maxEncodeBaseLen=19;
	  else maxEncodeBaseLen=len;
   }
	    
   public int getMaxEncodeBaseLen(){ 
	  return maxEncodeBaseLen;
   }
   
   public void setEncodeIntervalNum(int num){ 
	   encodeIntervalNum=num;
   }
	    
   public int getEncodeIntervalNum(){ 
	  return encodeIntervalNum;
   }
   
   public void setSplitStep(int step){ 
	   seqSplitStep=step;
   }
	    
   public int getSplitStep(){ 
	  return seqSplitStep;
   }
   
   public void filter(String inSeqFile, boolean isDeDup, int minSeqLen, String outDir){ 
	   
	   if(isDeDup) {
		  inSeqFile=removeDup(inSeqFile);
	   }
	   
	   if(minSeqLen>0) {
		  inSeqFile=SeqOperation.filterSeqFile(inSeqFile, seqSplitStep, minSeqLen, outDir);
	   }

   }
   
   public void filter(String inSeqFile, String inSeqFile2, boolean isDeDup,int minSeqLen, String outDir){ 
	   
	   String []inSeqFile12=new String[2];
	   if(isDeDup) {
		 inSeqFile12=removeDup(inSeqFile,inSeqFile2);
		 inSeqFile=inSeqFile12[0];
		 inSeqFile2=inSeqFile12[1];
	   }
	   
	   if(minSeqLen>0) {
		 inSeqFile12=SeqOperation.filterSeqFile(inSeqFile,inSeqFile2,seqSplitStep,minSeqLen,outDir);
	   }

   }

   public String removeDup(String inSeqFile){ 
	    
	    System.out.println("###### Cleaning duplicated single-end reads for:");	
	    System.out.println("Forward ["+inSeqFile+"]");
		String timeStamp=new SimpleDateFormat("yyyyMMdd_HHmmss").format(Calendar.getInstance().getTime());
		System.out.println("Start at time:"+timeStamp);
	
		System.out.println("Encoding reads......");
		List<SeqEncode> seqEncodeObjList=SeqOperation.getSeqNumEncode(inSeqFile,seqSplitStep,maxEncodeBaseLen);	
		timeStamp=new SimpleDateFormat("yyyyMMdd_HHmmss").format(Calendar.getInstance().getTime());
		System.out.println("Totally encoded reads: "+seqEncodeObjList.size()+"[time: "+timeStamp+"]");	
		
		System.out.println("Dispatching reads (pair-end)......");
		String outTmpDir=GEAT.getTmpDir()+"/"+timeStamp;
		FileOperation.newFolder(outTmpDir);
		List<String>dividedSeqFiles=SeqOperation.dispatchSeqByNumEncode(inSeqFile,seqSplitStep,
				seqEncodeObjList,encodeIntervalNum,outTmpDir);
		seqEncodeObjList=null;
		timeStamp=new SimpleDateFormat("yyyyMMdd_HHmmss").format(Calendar.getInstance().getTime());
		System.out.println("Reads dispatch done! "+"[time: "+timeStamp+"]");	
		
		List<SeqInfo> seqObjList=null;
	    List<SeqInfo> seqDupObjList;  
	    List<SeqInfo> seqCleanObjList;  
		int totalClean=0;
		int totalDup=0;
	    String seqFormat=FileOperation.getFileFormat(inSeqFile);
		String cleanSeqFile=inSeqFile.substring(0,inSeqFile.lastIndexOf("."))+".dupClean."+seqFormat;
		String dupSeqFile=inSeqFile.substring(0,inSeqFile.lastIndexOf("."))+".dup."+seqFormat;
		try {
			System.out.println("Cleaning duplicates...... ");	
			BufferedWriter writerClean=new BufferedWriter(new FileWriter(cleanSeqFile));
			BufferedWriter writerDup=new BufferedWriter(new FileWriter(dupSeqFile));
			int num=0;
			
			for(String seqFile:dividedSeqFiles) {		
			   num++;
			   seqObjList=SeqOperation.getSeqObj(seqFile);
			   String seqEncodeFile=seqFile.substring(0,seqFile.lastIndexOf("."))+".encode";
			   List<ArrayList <String>> encodeList=FileOperation.getMatrixFromFile(seqEncodeFile);
			   if(encodeList.size()==seqObjList.size()) {
				  int s=0;
				  while(s<seqObjList.size()) {
					 if(seqObjList.get(s).seqName.equals(encodeList.get(s).get(0))) {
					    seqObjList.get(s).seqNumEncode=Long.parseLong(encodeList.get(s).get(1));
					    seqObjList.get(s).seqNumRevEncode=Long.parseLong(encodeList.get(s).get(2));
					    s++;
					 }else {
					    System.out.println("Warning: seq name "+seqObjList.get(s).seqIdentifier
							   +" can not match each other["+encodeList.get(s).get(1)+"], just skip it!!!");
					    seqObjList.remove(s);
					 }
				  }
				  
				  setDupOfSeqObj(seqObjList);
				  
				  seqDupObjList= new ArrayList<SeqInfo>();  
				  seqCleanObjList= new ArrayList<SeqInfo>();  
				  for(SeqInfo seq: seqObjList) {
					 if(seq.isDup) seqDupObjList.add(seq);
					 else seqCleanObjList.add(seq);
				  }
				  
				  totalClean=totalClean+seqCleanObjList.size();	
				  totalDup=totalDup+seqDupObjList.size();
					///*
				  if(seqCleanObjList.size()>0) {
					if(SeqOperation.isFASTAFormat(cleanSeqFile)) {
					  SeqOperation.saveSeqObjAsFASTA(seqCleanObjList, writerClean);		
					}else if(SeqOperation.isFASTQFormat(cleanSeqFile)) {
					  SeqOperation.saveSeqObjAsFASTQ(seqCleanObjList, writerClean);  
					}
				  }
					
				  if(seqDupObjList.size()>0) {
					if(SeqOperation.isFASTAFormat(dupSeqFile)) {				  
					  SeqOperation.saveSeqObjAsFASTA(seqDupObjList, writerDup);
					}else if(SeqOperation.isFASTQFormat(dupSeqFile)) {
					  SeqOperation.saveSeqObjAsFASTQ(seqDupObjList, writerDup);  
					}	
				  }
					//*/		
				  seqCleanObjList=null;
				  seqDupObjList=null;				  
			   }else {
				  System.out.println("Warning: The size of NumEncode list is not equal to one of seq list, just skip it!!!");	
			   }
			   encodeList=null;
			   seqEncodeFile=null;
			   seqObjList=null;
			   timeStamp=new SimpleDateFormat("yyyyMMdd_HHmmss").format(Calendar.getInstance().getTime());
			   System.out.println(100.0f*(num*1.0f)/(dividedSeqFiles.size()*1.0f)+"% done! [time: "+timeStamp+"]");	
			}
			
			writerClean=null;
			writerDup=null;
			dividedSeqFiles=null;
			
			System.out.println("Reads after cleaning duplicates: "+totalClean);	
			System.out.println("Duplicated reads: "+totalDup);		
			System.out.println("Clean seq were saved as ["+cleanSeqFile+"]");	
			System.out.println("Duplicated seq were saved as ["+dupSeqFile+"]");
		    timeStamp=new SimpleDateFormat("yyyyMMdd_HHmmss").format(Calendar.getInstance().getTime());
			System.out.println("Duplicates clean done! [time: "+timeStamp+"]");	
		}catch(Exception e){
			System.out.println(e);
		}
		
	    FileOperation.delFolder(outTmpDir);
	    
	    return cleanSeqFile;

   }
   
   public String [] removeDup(String inSeqFile,String inSeqFile2){ 		
	    
	    System.out.println("###### Cleaning duplicated pair-end reads for:");
	    System.out.println("Forward ["+inSeqFile+"]");
	    System.out.println("Reverse ["+inSeqFile2+"]");
		String timeStamp=new SimpleDateFormat("yyyyMMdd_HHmmss").format(Calendar.getInstance().getTime());
		System.out.println("Start at time:"+timeStamp);
	
		System.out.println("Encoding reads......");
		List<SeqEncode> seqEncodeObjList=SeqOperation.getSeqNumEncode(inSeqFile,seqSplitStep,maxEncodeBaseLen);	
		timeStamp=new SimpleDateFormat("yyyyMMdd_HHmmss").format(Calendar.getInstance().getTime());
		System.out.println("Totally Encoded reads: "+seqEncodeObjList.size()+"[time: "+timeStamp+"]");	
	
		System.out.println("Dispatching reads(pair-end)......");
		String outTmpDir=GEAT.getTmpDir()+"/"+timeStamp;
		FileOperation.newFolder(outTmpDir);
		List<ArrayList<String>> dividedSeqFiles12=SeqOperation.dispatchSeqByNumEncode(inSeqFile,inSeqFile2,
				seqSplitStep,seqEncodeObjList,encodeIntervalNum,outTmpDir);
		seqEncodeObjList=null;

		timeStamp=new SimpleDateFormat("yyyyMMdd_HHmmss").format(Calendar.getInstance().getTime());
		System.out.println("Reads dispatch done! "+"[time: "+timeStamp+"]");	
		
		List<SeqInfo> seqObjList=null;
	    List<SeqInfo> seqDupObjList;  
	    List<SeqInfo> seqCleanObjList;  
		List<SeqInfo> seqObjList2=null;
	    List<SeqInfo> seqDupObjList2;  
	    List<SeqInfo> seqCleanObjList2; 
		int totalClean=0;
		int totalDup=0;
		
	    String seqFormat=FileOperation.getFileFormat(inSeqFile);
		String cleanSeqFile=inSeqFile.substring(0,inSeqFile.lastIndexOf("."))+".dupClean."+seqFormat;
		String dupSeqFile=inSeqFile.substring(0,inSeqFile.lastIndexOf("."))+".dup."+seqFormat;
		
		seqFormat=FileOperation.getFileFormat(inSeqFile2);
		String cleanSeqFile2=inSeqFile2.substring(0,inSeqFile2.lastIndexOf("."))+".dupClean."+seqFormat;
		String dupSeqFile2=inSeqFile2.substring(0,inSeqFile2.lastIndexOf("."))+".dup."+seqFormat;
		
		try {
			System.out.println("Cleaning duplicates...... ");	
			BufferedWriter writerClean=new BufferedWriter(new FileWriter(cleanSeqFile));
			BufferedWriter writerDup=new BufferedWriter(new FileWriter(dupSeqFile));
			BufferedWriter writerClean2=new BufferedWriter(new FileWriter(cleanSeqFile2));
			BufferedWriter writerDup2=new BufferedWriter(new FileWriter(dupSeqFile2));
			int num=0;			
			for(int d=0;d<dividedSeqFiles12.get(0).size();d++) {
			   num++;
			   String seqFile=dividedSeqFiles12.get(0).get(d);			  
			   seqObjList=SeqOperation.getSeqObj(seqFile);
			   String seqEncodeFile=seqFile.substring(0,seqFile.lastIndexOf("."))+".encode";
			   List<ArrayList <String>> encodeList=FileOperation.getMatrixFromFile(seqEncodeFile);
			   
			   String seqFile2=dividedSeqFiles12.get(1).get(d);
			   seqObjList2=SeqOperation.getSeqObj(seqFile2);
			  
			   //System.out.println("Checking pair-end seq......");
			   if(SeqOperation.checkPairEndSeq(seqObjList,seqObjList2)){			
			     if(encodeList.size()==seqObjList.size()){
				    int s=0;
				    while(s<seqObjList.size()) {
					  //forward & reverse
					  if(seqObjList.get(s).seqName.equals(encodeList.get(s).get(0))) {						  
					     seqObjList.get(s).seqNumEncode=Long.parseLong(encodeList.get(s).get(1));
					     seqObjList.get(s).seqNumRevEncode=Long.parseLong(encodeList.get(s).get(2));				    
					     s++;
					  }else {
					     System.out.println("Warning: forward seq name "+seqObjList.get(s).seqIdentifier
							   +" can not match each other["+encodeList.get(s).get(0)+"], just skip it!!!");
					     seqObjList.remove(s);
					     seqObjList2.remove(s);
					  }
				    }
				  
				    setDupOfSeqObj(seqObjList,seqObjList2);				  

				    seqDupObjList= new ArrayList<SeqInfo>();  
				    seqCleanObjList= new ArrayList<SeqInfo>();  
				    seqDupObjList2= new ArrayList<SeqInfo>();  
				    seqCleanObjList2= new ArrayList<SeqInfo>(); 
				   
				    for(int i=0;i<seqObjList.size();i++) {
					  if(seqObjList.get(i).isDup) {
					    seqDupObjList.add(seqObjList.get(i));
						seqDupObjList2.add(seqObjList2.get(i));
					  }else {
						seqCleanObjList.add(seqObjList.get(i));
						seqCleanObjList2.add(seqObjList2.get(i));
					  }
				    }
				  
				    totalClean=totalClean+seqCleanObjList.size();	
				    totalDup=totalDup+seqDupObjList.size();
				   
				    ///*
				    //save seq
				    //forward
				    if(seqCleanObjList.size()>0) {
					  if(SeqOperation.isFASTAFormat(cleanSeqFile)) {
					    SeqOperation.saveSeqObjAsFASTA(seqCleanObjList, writerClean);		
					  }else if(SeqOperation.isFASTQFormat(cleanSeqFile)) {
					    SeqOperation.saveSeqObjAsFASTQ(seqCleanObjList, writerClean);  
					  }
				    }
					
				    if(seqDupObjList.size()>0) {
					  if(SeqOperation.isFASTAFormat(dupSeqFile)) {				  
					    SeqOperation.saveSeqObjAsFASTA(seqDupObjList, writerDup);
					  }else if(SeqOperation.isFASTQFormat(dupSeqFile)) {
					    SeqOperation.saveSeqObjAsFASTQ(seqDupObjList, writerDup);  
					  }	
				    }
				   
				    //reverse
				    if(seqCleanObjList2.size()>0) {
					  if(SeqOperation.isFASTAFormat(cleanSeqFile2)) {
					    SeqOperation.saveSeqObjAsFASTA(seqCleanObjList2, writerClean2);		
					  }else if(SeqOperation.isFASTQFormat(cleanSeqFile2)) {
					    SeqOperation.saveSeqObjAsFASTQ(seqCleanObjList2, writerClean2);  
					  }
				    }
					
				    if(seqDupObjList2.size()>0) {
					  if(SeqOperation.isFASTAFormat(dupSeqFile2)) {				  
					    SeqOperation.saveSeqObjAsFASTA(seqDupObjList2, writerDup2);
					  }else if(SeqOperation.isFASTQFormat(dupSeqFile2)) {
					    SeqOperation.saveSeqObjAsFASTQ(seqDupObjList2, writerDup2);  
					  }	
				    }
					//*/		
				    seqCleanObjList=null;
				    seqDupObjList=null;
				    seqCleanObjList2=null;
				    seqDupObjList2=null;				  
			     }else {
				    System.out.println("Warning: The encode list size is not equal to seq list, just skip it!!!");	
			     }
			     encodeList=null;
			     seqEncodeFile=null;
			     seqObjList=null;
			     seqObjList2=null;
			     timeStamp=new SimpleDateFormat("yyyyMMdd_HHmmss").format(Calendar.getInstance().getTime());
			     System.out.println(100.0f*(num*1.0f)/(dividedSeqFiles12.get(0).size()*1.0f)+"% done! [time: "+timeStamp+"]");	
			  }else { 
				 System.err.println("Warning: Dispatch No."+num+" of pair-end reads are not paired each other. Just skip it!!");
			  }			 
			}
			
			writerClean=null;
			writerDup=null;
			writerClean2=null;
			writerDup2=null;
			dividedSeqFiles12=null;
			
			System.out.println("Reads after cleaning duplicates: "+totalClean);	
			System.out.println("Duplicated reads: "+totalDup);		
			System.out.println("Clean seq were saved as ["+cleanSeqFile+"]");	
			System.out.println("Duplicated seq were saved as ["+dupSeqFile+"]");
		    timeStamp=new SimpleDateFormat("yyyyMMdd_HHmmss").format(Calendar.getInstance().getTime());
			System.out.println("Duplicates clean done! [time: "+timeStamp+"]");	
		}catch(Exception e){
			System.out.println(e);
		}
		
	    FileOperation.delFolder(outTmpDir);
	    
	    String []inSeqFile12=new String[2];
	    inSeqFile12[0]=cleanSeqFile;
	    inSeqFile12[1]=cleanSeqFile2;
	    
	    return inSeqFile12;

   }
 
   void setDupOfSeqObj(List<SeqInfo> seqObjList) {
	    
	    //String timeStamp;
	    //System.out.println("Sorting reads......");
		boolean desc=false;
		Collections.sort(seqObjList, new SeqInfo.CompSeqEncode(desc));
		//System.out.println("Sorted reads: "+seqObjList.size());
		//timeStamp=new SimpleDateFormat("yyyyMMdd_HHmmss").format(Calendar.getInstance().getTime());
		//System.out.println("time: "+timeStamp);
	   
		SeqInfo seqObj;
		SeqInfo seqObj_r0 = null;
		SeqInfo currSeqObj;
		SeqInfo currSeqObj2;  
	    List<SeqInfo> querySeqObjList;  
		
	    int i=0;  
	    int r=0;
	    //int m=1;
		for(int s=0;s<seqObjList.size();s++){
		  currSeqObj=seqObjList.get(s);	
		  if(!currSeqObj.isDup && currSeqObj.dupNum<0){
			  /*
			 if((s+1)/500000==m) {
				timeStamp=new SimpleDateFormat("yyyyMMdd_HHmmss").format(Calendar.getInstance().getTime());
				System.out.println((s+1)+" time: "+timeStamp);
				m++;
			 }	
			 */
			 querySeqObjList=SeqOperation.findSeqObjOfSameEncode(seqObjList,currSeqObj);
			 if(querySeqObjList.size()==1) {
				currSeqObj.isDup=false;
				currSeqObj.dupNum=0; 
			 }else if(querySeqObjList.size()>1){
			    currSeqObj2=currSeqObj;
			    while(currSeqObj2!=null){		 
			      currSeqObj2=null;
				  r=0;
				  i=0;
				  while(i<querySeqObjList.size()) {	
				    seqObj=querySeqObjList.get(i);		 
				    if(seqObj.seq.equalsIgnoreCase(currSeqObj.seq)){				
					  if(r==0) { 			
					    seqObj_r0=seqObj;
					    seqObj_r0.isDup=false;
					  }if(r>0) {				    
					    seqObj.isDup=true;
					    seqObj.dupNum=r;	
					  }
					  seqObj_r0.dupNum=r;
					  querySeqObjList.remove(seqObj);
					  r++;
				    }else{
					  currSeqObj2=seqObj;	
					  i++;		
				    } 		    
				  }// while query
				  if(currSeqObj2!=null) currSeqObj=currSeqObj2;
			    }
			 }
			 querySeqObjList=null;	  
		  }	  
		}
		
	    Collections.sort(seqObjList, new SeqInfo.CompSeqID(desc));	
   }
   
   void setDupOfSeqObj(List<SeqInfo> seqObjList, List<SeqInfo> seqObjList2) {
		
	    //System.out.println("Sorting reads......");
		boolean desc=false;
		Collections.sort(seqObjList, new SeqInfo.CompSeqEncode(desc));		
		//String timeStamp=new SimpleDateFormat("yyyyMMdd_HHmmss").format(Calendar.getInstance().getTime());
		//System.out.println("Sorted reads: "+seqObjList.size()+"[time: "+timeStamp+"]");		
		
		SeqInfo seqObj;
		SeqInfo seqObj_r0 = null;
		SeqInfo currSeqObj;
		SeqInfo currSeqObj2;  
	    List<SeqInfo> querySeqObjList; 	
	    int i=0;  
	    int r=0;
	    //int m=1;
		for(int s=0;s<seqObjList.size();s++){
		  currSeqObj=seqObjList.get(s);	
		  if(!currSeqObj.isDup && currSeqObj.dupNum<0){
			  /*
			 if((s+1)/500000==m) {
				timeStamp=new SimpleDateFormat("yyyyMMdd_HHmmss").format(Calendar.getInstance().getTime());
				System.out.println((s+1)+" time: "+timeStamp);
				m++;
			 }	
			 */
			 querySeqObjList=SeqOperation.findSeqObjOfSameEncode(seqObjList,currSeqObj);
			 if(querySeqObjList.size()==1) {
				currSeqObj.isDup=false;
				currSeqObj.dupNum=0; 
			 }else if(querySeqObjList.size()>1){
			    currSeqObj2=currSeqObj;
			    while(currSeqObj2!=null){		 
			      currSeqObj2=null;
				  r=0;
				  i=0;
				  while(i<querySeqObjList.size()){	
				    seqObj=querySeqObjList.get(i);		 
				    if(seqObj.seq.equalsIgnoreCase(currSeqObj.seq) &&
				       seqObjList2.get(seqObj.seqID).seq.equalsIgnoreCase(seqObjList2.get(currSeqObj.seqID).seq)){				
					   if(r==0){ 				
						 seqObj_r0=seqObj;
						 seqObj_r0.isDup=false;
					   }if(r>0){				     
					     seqObj.isDup=true;
					     seqObj.dupNum=r;			
					   }
					   seqObj_r0.dupNum=r;
					   querySeqObjList.remove(seqObj);
					   r++;
				    }else{
					   currSeqObj2=seqObj;	
					   i++;		
				    } 		    
				  }// while query
				  if(currSeqObj2!=null) currSeqObj=currSeqObj2;
			    }
			 }
			 querySeqObjList=null;	  
		  }	  
		}
		
	    Collections.sort(seqObjList, new SeqInfo.CompSeqID(desc));	
	    Collections.sort(seqObjList2, new SeqInfo.CompSeqID(desc));	
	   
   }
 
    
   public static int removeSeqAndItsDup(String seqFile,String seqNameFile,int seqNameColIdx){
	
	    int num=0;
		List<SeqInfo> subSeqObjList=SeqOperation.extractSubSeq(seqFile,seqNameFile,seqNameColIdx);
		
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
		String format=FileOperation.getFileFormat(seqFile);
		String outSeqFile=seqNameFile+"Dup."+format;
		SeqOperation.saveSeqObj(dupSeqObjList, outSeqFile);
		dupSeqObjList=null;
		 
		List<SeqInfo> filterSeqObjList=SeqOperation.excludeSeq(seqObjList,excludedSeqNameList);
		outSeqFile=seqFile.substring(0,seqFile.lastIndexOf("."))+".clean."+format;
		SeqOperation.saveSeqObj(filterSeqObjList, outSeqFile);
		 
		num=filterSeqObjList.size();
		 
		filterSeqObjList=null;
		 
		return num;
   }
 

}