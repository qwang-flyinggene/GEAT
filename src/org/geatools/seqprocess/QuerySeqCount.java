package org.geatools.seqprocess;
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

import java.io.File;
import java.text.SimpleDateFormat;
import java.util.ArrayList; 
import java.util.Calendar;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;

import org.geatools.data.structure.SeqCompoAlignInfo;
import org.geatools.data.structure.SeqCompoRecognizer;
import org.geatools.data.structure.SeqCompoFeatures;
import org.geatools.data.structure.SeqInfo;
import org.geatools.data.structure.SeqRocket;
import org.geatools.operation.FileOperation;
import org.geatools.operation.SeqOperation;

public class  QuerySeqCount extends SeqRocketRecognition{  
  
  int maxEncodeBaseLen=19;
  int localExactEncodeLen=3;
  //float nonexactRate=0.7f;
  boolean doNonExact=true;
  
  int minQueryLen=2;  
  double queryMinAlignRatio=0.9d;
  double queryMaxMismatchRatio=0.1d;
  double queryMaxGapRatio=0.1d;
  String querySeqName="QuerySeq";
  String querySeqColor="red"; 

  
  final HashMap<Character, Integer> base2numMap =new HashMap<Character, Integer>();  
  
  public QuerySeqCount(){ 

  }
  public QuerySeqCount(String homeDir,String dataDir,String tmpDir){ 
	  setHomeDir(homeDir);
	  setDataDir(dataDir);
	  setTmpDir(tmpDir);	 
  }
 
  public void setDoNonexact(boolean action){
	  doNonExact=action;
  }
  
  public void setMaxEncodeBaseLen(int len){ 
	  maxEncodeBaseLen=len;
  }
    
  public int getMaxEncodeBaseLen(){ 
	  return maxEncodeBaseLen;
  }

  public void setLocalExactEncodeLen(int len){
	  localExactEncodeLen=len;
  }
  
  void configRecognizer( List<SeqRocket> rockets){
    
	SeqCompoRecognizer query;
	for(int i=0;i<rockets.size();i++){	  
	  query=null;	  
      SeqCompoRecognizer recognizer;
	  String seqType="";
	  for(int k=0;k<rockets.get(i).seqRecognizers.size();k++){
	    recognizer=rockets.get(i).seqRecognizers.get(k);	
	    seqType=rockets.get(i).seqCompoFeatures.compoNames.get(recognizer.index);
	    if(seqType.equals(querySeqName)) query=recognizer;
	  }
	  recognizer=null;
	  
	  if(query!=null && query.rawSeq!=null && !query.rawSeq.equals("")){

		  query.seq=query.rawSeq;	
		  query.seqLength=query.seq.length();
		  query.territoryLen=query.seqLength+(int) Math.ceil(query.seqLength*territoryLeftExtendRatio);	
		  query.exactMaxStart=query.territoryLen-query.seqLength+1;
		  query.minAlignRatio=queryMinAlignRatio;
		  query.maxMismatchRatio=queryMaxMismatchRatio;
		  query.maxGapRatio=queryMaxGapRatio;			
		  query.minAlignLen=(int) Math.ceil(query.seqLength*query.minAlignRatio);
		  query.maxMismatchNum=(int) Math.ceil(query.seqLength*query.maxMismatchRatio);
		  query.maxGapNum=(int) Math.ceil(query.seqLength*query.maxGapRatio);	
		  if(query.minAlignLen<=9){	  
			  query.blastWordSize=4;
			  query.blastTask="blastn-short";
		  }else if(query.minAlignLen<=16){
			  query.blastWordSize=7;
			  query.blastTask="blastn-short";
		  }else if(query.minAlignLen<=26){
			  query.blastWordSize=11;
			  query.blastTask="blastn-short";
		  }else if(query.minAlignLen<=40){
			  query.blastWordSize=16;
			  query.blastTask="blastn-short";
	      }else if(query.minAlignLen<=64){
	    	  query.blastWordSize=24;
	    	  query.blastTask="blastn";
	      }else{
	    	  query.blastWordSize=32;
	    	  query.blastTask="megablast";
		  }
		  
		  if(query.seqLength<4){
			  query.maxQStart=1;
			  query.maxSStart=1;
		  }else if(query.seqLength<minQueryLen){
			  query.maxQStart=query.seqLength/2;
			  query.maxSStart=query.seqLength/2;
          }else{			
        	  query.maxQStart=query.seqLength-query.minAlignLen+1;
        	  query.maxSStart=query.territoryLen-query.minAlignLen+1;
		  }	   
		  query.leftSubForBLAST=false;
		  query.leftShiftForNext=false;
		  //query.saveRecognizedSeq=false;  //ooooooooooooooo
		  query.saveAsFASTA=false; //ooooooooooooooo	
		  query.saveAsFASTQ=false; //ooooooooooooooo	
		  query.leftTrimSave=false;
		  query.rightTrimSave=false;
		  query.exactAlignedSeqFile=tmpDir+File.separator+query.seqName+"_ExactAlignedSeq";
		  query.seqFASTAFile=tmpDir+File.separator+query.seqName+".fna";		  

          rockets.get(i).seqRecognizers.set(query.index,query);	

      }else if(query.rawSeq==null || query.rawSeq.equals("")){
    	  rockets.get(i).seqRecognizers.set(query.index,null);	
      }
	 	  
	  query=null;
	 
	}	
   
  } 
  
  /*
  List<SeqCompoAlignInfo> getNonExactSeq(List<SeqCompoAlignInfo> seqObj, 
		  List<SeqRocket> rockets, String inSeqFileFormat, String seqType){
		
		 List<SeqCompoAlignInfo> nonExactMatchedSeqList=new ArrayList<SeqCompoAlignInfo>();
	
		 String seqLine;	
		 boolean isMatch=false;	
		 int seqStartIdx=0;		
		 SeqCompoRecognizer querySeq=null;
		 for(int i=0;i<seqObj.size();i++){ 
				seqLine=seqObj.get(i).seq;				
				isMatch=false;	
				for(int j=0;j<rockets.size();j++){
				
					querySeq=rockets.get(j).seqRecognizer.get(
							 rockets.get(j).seqTypeInfo.seqTypeName.indexOf(seqType));
					seqStartIdx=seqLine.indexOf(querySeq.seq);			 
					if(seqStartIdx>=0 && seqStartIdx<=querySeq.maxExactStart-1){
					  isMatch=true;
					  rockets.get(j).exactSeqCounts++;
					  
					  break;
					} 
				}
				
				if(!isMatch){						 
			      nonExactMatchedSeqList.add(seqObj.get(i));				  
				}
		 }

		 seqLine=null;
		
		 return nonExactMatchedSeqList;

  }
  */
  
  void setExactAlignInfo(List<SeqCompoAlignInfo> sortedSeqList,List<SeqRocket> rockets,String seqCompName){		
		 
	     if(sortedSeqList==null || rockets==null || rockets.size()==0) return;		
		 
		 for(SeqRocket seqRocket:rockets){				 
			 int leftSStartIndex=0;
		     int leftSStart=0;
		     int leftSEnd=0;				
			 int trimSEndIndex=-1;
			 int trimLeftShift=0;
			 int trimSEnd=0;			
			 for(int i=0;i<seqRocket.seqRecognizers.size();i++){
			   if(seqRocket.seqRecognizers.get(i).done 
					  && seqRocket.seqRecognizers.get(i).leftShiftForNext ){		 
					    
				   trimSEndIndex=seqRocket.seqRecognizers.get(i).index;
				   trimLeftShift=seqRocket.seqRecognizers.get(i).trimLeftShift;					  
			   }
			 }					
			 int alignSStartIndex=seqRocket.seqCompoFeatures.compoNames.indexOf(seqCompName);		  
			 int alignSEndIndex=alignSStartIndex;

			 SeqCompoRecognizer querySeq=seqRocket.seqRecognizers.get(
					 seqRocket.seqCompoFeatures.compoNames.indexOf(seqCompName));	
		     int leftMaxStartIndex=querySeq.exactMaxStart-1;
		     if (leftMaxStartIndex<0) leftMaxStartIndex=0;	
			 int seqIdx = Collections.binarySearch(
					 sortedSeqList,
					 querySeq, 
					 new SeqCompoRecognizer.CompSeqEncode(false)
			 );
			 
			 if(seqIdx>=0){
			   if(trimSEndIndex<0){
				 //go backward......
				 for(int j=seqIdx;j>=0;j--){				
				   ////leftSStartIndex=sortedSeqList.get(j).seq.indexOf(querySeq.seq); 
				   ////if(leftSStartIndex>=0 && leftSStartIndex<=leftMaxStartIndex){
				   if(sortedSeqList.get(j).seq.equalsIgnoreCase(querySeq.seq)){
					   //if( sortedSeqList.get(j).seqAlignSStart.get(alignSStartIndex)<0){
	  				  seqRocket.exactSeqCount++;					 
					  ////leftSStart=leftSStartIndex+1;
					  ////leftSEnd=leftSStartIndex+querySeq.seq.length();
	  				  leftSStart=1;
					  leftSEnd=querySeq.seq.length();
					  sortedSeqList.get(j).seqAlignSStart.set(alignSStartIndex,leftSStart);
					  sortedSeqList.get(j).seqAlignSEnd.set(alignSEndIndex,leftSEnd);
					   //}else{
						   //System.out.println("oooooo: "+querySeq.seqName+" "+querySeq.seq+" "+querySeq.seqNumEncode);
					   //}
				   }
				   if(querySeq.seqNumEncode!=sortedSeqList.get(j).seqNumEncode) break;				   
				 }
				 //go forward......
			     for(int j=seqIdx+1;j<sortedSeqList.size();j++){				
				   ////leftSStartIndex=sortedSeqList.get(j).seq.indexOf(querySeq.seq); 
				   ////if(leftSStartIndex>=0 && leftSStartIndex<=leftMaxStartIndex){
				   if(sortedSeqList.get(j).seq.equalsIgnoreCase(querySeq.seq)){
					   //if( sortedSeqList.get(j).seqAlignSStart.get(alignSStartIndex)<0){
					  seqRocket.exactSeqCount++;					 
					  ////leftSStart=leftSStartIndex+1;
					  ////leftSEnd=leftSStartIndex+querySeq.seq.length();
	  				  leftSStart=1;
					  leftSEnd=querySeq.seq.length();
					  sortedSeqList.get(j).seqAlignSStart.set(alignSStartIndex,leftSStart);
					  sortedSeqList.get(j).seqAlignSEnd.set(alignSEndIndex,leftSEnd);
					   //}else{
						   //System.out.println("oooooo: "+querySeq.seqName+" "+querySeq.seq+" "+querySeq.seqNumEncode);
					   //}
				   }
				   if(querySeq.seqNumEncode!=sortedSeqList.get(j).seqNumEncode) break;				   
			     }
			   }else{
				 //go backward......
				 for(int j=seqIdx;j>=0;j--){	
				   leftSStartIndex=sortedSeqList.get(j).seq.indexOf(querySeq.seq);		
				   trimSEnd=sortedSeqList.get(j).seqAlignSEnd.get(trimSEndIndex)-trimLeftShift;		   	
				   //if(leftSStartIndex>=(trimSEnd-1) && leftSStartIndex<=(trimSEnd+leftMaxStartIndex)){
				   if(sortedSeqList.get(j).seq.length()==(trimSEnd+querySeq.seq.length())
						   && leftSStartIndex>=(trimSEnd-1)){
				      seqRocket.exactSeqCount++;
					  leftSStart=leftSStartIndex+1;
					  leftSEnd=leftSStartIndex+querySeq.seq.length();		
					  sortedSeqList.get(j).seqAlignSStart.set(alignSStartIndex,leftSStart);
					  sortedSeqList.get(j).seqAlignSEnd.set(alignSEndIndex,leftSEnd);	
				   }
				   if(querySeq.seqNumEncode!=sortedSeqList.get(j).seqNumEncode) break;
				 }
				 //go forward......
				 for(int j=seqIdx+1;j<sortedSeqList.size();j++){	
				   leftSStartIndex=sortedSeqList.get(j).seq.indexOf(querySeq.seq);		
				   trimSEnd=sortedSeqList.get(j).seqAlignSEnd.get(trimSEndIndex)-trimLeftShift;		   	
				   //if(leftSStartIndex>=(trimSEnd-1) && leftSStartIndex<=(trimSEnd+leftMaxStartIndex)){
				   if(sortedSeqList.get(j).seq.length()==(trimSEnd+querySeq.seq.length())
						   && leftSStartIndex>=(trimSEnd-1)){
					  seqRocket.exactSeqCount++;
					  leftSStart=leftSStartIndex+1;
					  leftSEnd=leftSStartIndex+querySeq.seq.length();		
					  sortedSeqList.get(j).seqAlignSStart.set(alignSStartIndex,leftSStart);
					  sortedSeqList.get(j).seqAlignSEnd.set(alignSEndIndex,leftSEnd);	
				   }
				   if(querySeq.seqNumEncode!=sortedSeqList.get(j).seqNumEncode) break;
				 }

			   }
		     }else{
               System.out.println("NoMatch: "+querySeq.seqName+" "+querySeq.seq+" "+querySeq.seqNumEncode);
		     }
			 querySeq=null;
		 }
  }
  
  List<SeqCompoAlignInfo> findQuerySeq(List<SeqCompoAlignInfo> sortedSeqList,SeqCompoRecognizer querySeq){
		 
	     List<SeqCompoAlignInfo> querySeqList= new ArrayList<SeqCompoAlignInfo>();
	     int seqIdx = Collections.binarySearch(
	    		 sortedSeqList,
				 querySeq, 
				 new SeqCompoRecognizer.CompSeqEncode(false)
		 );

		 if(seqIdx>=0){		  
			 //go backward......
			 for(int j=seqIdx;j>=0;j--){			   
			   if(querySeq.seqNumEncode!=sortedSeqList.get(j).seqNumEncode) break;
			  
			   if(querySeq.seqNumEncode==sortedSeqList.get(j).seqNumEncode
					&& querySeq.seqNumRevEncode==sortedSeqList.get(j).seqNumRevEncode)
			   querySeqList.add(sortedSeqList.get(j));
			 }
			 //go forward......
		     for(int j=seqIdx+1;j<sortedSeqList.size();j++){		     
		       if(querySeq.seqNumEncode!=sortedSeqList.get(j).seqNumEncode) break;
			   
			   if(querySeq.seqNumEncode==sortedSeqList.get(j).seqNumEncode
					&& querySeq.seqNumRevEncode==sortedSeqList.get(j).seqNumRevEncode)
		       querySeqList.add(sortedSeqList.get(j));
		     }
		 }
		 
		 return querySeqList;
  }
 
 
  public void launchQueryCount(List<String>expSeqFiles,String querySeqFile,
			 String outDir, String outTag, boolean doSplitSeq){	    
	  
	  List<List<SeqRocket>> expList=new ArrayList<List<SeqRocket>>();
	  // forward seq(Single-end)
	  for(String inSeqFile:expSeqFiles){
		 if(new File(inSeqFile).exists()){
		   launchQueryCount(inSeqFile, querySeqFile,outDir,null,true);			
		   expList.add(getSeqRockets());	
		 }else{
		   System.out.println("Warning: ["+inSeqFile+"] doesn't exist, we skipped it!!!");
		 }
	  }
	  
	  if(expList.size()==0) return;
	  
	  if(outDir==null){
		 outDir=querySeqFile.substring(0, querySeqFile.lastIndexOf(File.separator))+File.separator+"RecognizedSeq";	
	  }
	  FileOperation.newFolder(outDir);

	  //Save query seq count ................................
	  if(outTag!=null) outTag=outTag+"_";
	  else outTag="";
	  
	  String outFile=outTag+querySeqFile.substring(
					   querySeqFile.lastIndexOf(File.separator)+1,querySeqFile.lastIndexOf(".")
				    )+".Counts.txt";
      
	  outFile=outDir+File.separator+outFile;
	  List<ArrayList<String>> resOut=new ArrayList<ArrayList<String>>();
	  ArrayList<String> perRes;
  	  for(int i=0;i<expList.get(0).size();i++){	
		 perRes=new ArrayList<String>();
		 perRes.add(expList.get(0).get(i).rocketName);
		 for(List<SeqRocket> exp:expList){
		   perRes.add(Integer.toString(exp.get(i).recognizedSeqCount));	
		 }
		 resOut.add(perRes);
		 perRes=null;		 
	  }// for Rocket
  	  FileOperation.saveMatrixList(resOut, outFile);
      resOut=null;
      
      //Save query seq RPM ................................
	  outFile=outTag+querySeqFile.substring(
				querySeqFile.lastIndexOf(File.separator)+1,querySeqFile.lastIndexOf(".")
			  )+".RPM.txt";
      
	  outFile=outDir+File.separator+outFile;
	  resOut=new ArrayList<ArrayList<String>>();	 
  	  for(int i=0;i<expList.get(0).size();i++){	
		 perRes=new ArrayList<String>();
		 perRes.add(expList.get(0).get(i).rocketName);
		 for(List<SeqRocket> exp:expList){
		   perRes.add(Double.toString(exp.get(i).recognizedSeqRPM));	
		 }
		 resOut.add(perRes);
		 perRes=null;		 
	  }// for Rocket
  	  FileOperation.saveMatrixList(resOut, outFile);
      resOut=null;
  }
  
  public List<SeqRocket> launchQueryCount(String inSeqFile,String querySeqFile,
		 String outDir, String outFile, boolean saveRes){	    
		
	    setSeqRockets(null);    

		// forward seq(Single-end)
		List<SeqRocket> querySeqRockets =buildQueryRockets(querySeqFile);		
		if(querySeqRockets==null || querySeqRockets.size()==0) return querySeqRockets;
		
		tmpFiles=new  ArrayList<String>();		
		List<SeqCompoAlignInfo> seqObjList=null;
		System.out.println("Counting for ["+inSeqFile+"]............");
		//Setting and Checkpoint for forward Seq................
		System.out.println("Total reads: "+SeqOperation.getSeqNum(inSeqFile));
		
		System.out.println("Checking reads......");
		seqObjList=checkSeq(inSeqFile);
		if(seqObjList==null || seqObjList.size()==0){
		   System.out.println("Warning: Empty sequence for ["+inSeqFile+"]");
		   return querySeqRockets;
		}
		System.out.println("Checked reads: "+seqObjList.size());
		
		System.out.println("Encoding reads......");
		setSeqNumEncode(seqObjList,maxEncodeBaseLen);		
		System.out.println("Sorting reads......");
		boolean desc=false;
		Collections.sort(seqObjList, new SeqCompoAlignInfo.CompSeqEncode(desc));
        
		//initialize seq align array
		initSeqAlignArray(seqObjList,querySeqRockets.get(0).seqCompoFeatures.compoNames.size());

		//Recognize and set exact query seq align information including counts,alignment position...............................................	        
		System.out.println("Identifying exact match......");	
		setExactAlignInfo(seqObjList,querySeqRockets,querySeqName);		
		
		//Recognize and set non-exact query seq ...............................................
		if(doNonExact){
		  System.out.println("Identifying non-exact alignment......");
		  //Get non-exact query seq ...............................................
		  int compoIdx=querySeqRockets.get(0).seqCompoFeatures.compoNames.indexOf(querySeqName);
		  List<SeqCompoAlignInfo> nonExactSeqList=getNoRecognizedSeq(seqObjList,compoIdx);		
		  seqObjList=null;				
			
		  //Set non-exact query seq align info...............................................
		  setNonExactAlignInfo(querySeqRockets,querySeqName,nonExactSeqList,localExactEncodeLen);			
	      nonExactSeqList=getNoRecognizedSeq(nonExactSeqList,compoIdx);	
	      System.out.println("o>o "+nonExactSeqList.size());
		  /*
		  nonExactSeqList=getNoRecognizedSeq(nonExactSeqList,compoIdx);		  
		  if(nonExactSeqList.size()>0){
			seqEncodeLen=3;
			setNonExactAlignInfo(querySeqRockets,querySeqName,nonExactSeqList,seqEncodeLen);
		  }
		  */
		  nonExactSeqList=null;  		  
		}
		
		//Set query seq count................................
		int totalQueryCount=0;
		for(SeqRocket seqRocket:querySeqRockets){
		  seqRocket.recognizedSeqCount=seqRocket.exactSeqCount+seqRocket.nonExactSeqCount;
		  totalQueryCount+=seqRocket.recognizedSeqCount;
		}
		for(SeqRocket seqRocket:querySeqRockets){	
		  seqRocket.recognizedSeqRPM=(Math.pow(10,6)*seqRocket.recognizedSeqCount)/(1.0d*totalQueryCount)+1;
		}// for Rocket	
		
		setSeqRockets(querySeqRockets);
		
		//System.out.println("Delete temporary files");
	    for(String tmpFile: tmpFiles){
	      FileOperation.delFile(tmpFile);
		}
	    	    
	    if(saveRes){
		   if(outDir==null){
				outDir=inSeqFile.substring(0, inSeqFile.lastIndexOf(File.separator))+File.separator+"RecognizedSeq";	
		   }
		   FileOperation.newFolder(outDir);

		   if(outFile==null)
				outFile=inSeqFile.substring(
						   inSeqFile.lastIndexOf(File.separator)+1,inSeqFile.lastIndexOf(".")
					    )+"_"
						+querySeqFile.substring(
						   querySeqFile.lastIndexOf(File.separator)+1,querySeqFile.lastIndexOf(".")
					    )+".Counts.txt";
           
		   outFile=outDir+File.separator+outFile;
		   saveSeqCountInfo(querySeqRockets,outFile);		   
	    }
	    
	    System.out.println("Done for ["+inSeqFile+"]");
	    
	    return querySeqRockets;
	   
  }
  
  List<SeqRocket> splitLaunchQueryCount(String inSeqFile,int splitStep,
			 String querySeqFile, String splitedSeqOut, String outDir,String outFile){
		 
		 List<String> splitedSeqFiles=new ArrayList<String>();
		//Check seq format, and then split seq into multiple subfiles................
		 System.out.println("Total Seq Num: "+SeqOperation.getSeqNum(inSeqFile));
		 splitedSeqFiles=SeqOperation.splitSeqFile(inSeqFile, splitStep, splitedSeqOut);	 
		 
		 List<SeqRocket> seqRockets=splitLaunchQueryCount(splitedSeqFiles,querySeqFile,outDir,outFile);
		 
		 return seqRockets;
  }
	 
  List<SeqRocket> splitLaunchQueryCount(List<String> splitedSeqFiles, String querySeqFile,
		  String outDir,String outFile){

		 List<SeqRocket> seqRockets=new ArrayList<SeqRocket>();
		 
		 List<ArrayList<SeqRocket>> splitedRockets=new ArrayList<ArrayList<SeqRocket>>();
		 ArrayList<SeqRocket> rockets;
		 
		 String splitOutDir = null;
		 int s=1;
		 for(String file:splitedSeqFiles){
		    try {
		      splitOutDir=file.substring(0,file.lastIndexOf(File.separator));
		      splitOutDir=splitOutDir+File.separator+file.substring(file.lastIndexOf(File.separator)+1,file.lastIndexOf("."));
			  FileOperation.newFolder(splitOutDir);
			   
			  System.out.println("......Recognizing sequence for split "+s+" ......");

			  rockets=(ArrayList<SeqRocket>) launchQueryCount(file,querySeqFile,splitOutDir,outFile,false);		  
		      if(rockets!=null && rockets.size()>0) splitedRockets.add(rockets);
			
			  s++;
			}catch (Exception e) {
			// TODO Auto-generated catch block
			   e.printStackTrace();
			}
		 }
		   
		 ///*
		 if(splitedSeqFiles.size()>0){
			 String file=splitedSeqFiles.get(0);
			 if(outDir==null) outDir=file.substring(0,file.lastIndexOf(File.separator))+File.separator+"combined";
			 FileOperation.newFolder(outDir);
			 System.out.println("......Combine splited sequences......");
			 int seqCounts=0;		
	         List<ArrayList<String>> res=new ArrayList<ArrayList<String>>();
			 ArrayList<String> per;
		     for(int r=0;r<splitedRockets.get(0).size();r++){
		    	 SeqRocket rocket=new SeqRocket();
		    	 rocket=splitedRockets.get(0).get(r);
				 per=new ArrayList<String>();
				 per.add(rocket.rocketName);
				 seqCounts=0;
				 for(int split=0;split<splitedRockets.size();split++){
					seqCounts=seqCounts+splitedRockets.get(split).get(r).recognizedSeqCount;
				 }
				 rocket.recognizedSeqCount=seqCounts;
				 per.add(Integer.toString(seqCounts));
				 res.add(per);
				 per=null;
				 
				 seqRockets.add(rocket);
				 rocket=null;
		     }
			 if(outFile==null)
			    outFile=querySeqFile.substring(
				  querySeqFile.lastIndexOf(File.separator)+1,querySeqFile.lastIndexOf("."))+".counts.txt";
			 
			 outFile=outDir+File.separator+outFile;
			 FileOperation.saveMatrixList(res, outFile);
		 }
		 //*/
		 splitedRockets=null;
		 
		 setSeqRockets(seqRockets);	 
		 System.out.println("......Seq Counts Done......");
		 
	     return seqRockets;
  }

  void saveSeqCountInfo(List<SeqRocket>querySeqRockets,String outFile){
	    //Save query seq count info................................
		List<ArrayList<String>> resOut=new ArrayList<ArrayList<String>>();
		ArrayList<String> perRes;
		for(SeqRocket seqRocket:querySeqRockets){	
			  perRes=new ArrayList<String>();
			  perRes.add(seqRocket.rocketName);
			  perRes.add(Integer.toString(seqRocket.recognizedSeqCount));
			  perRes.add(Double.toString(seqRocket.recognizedSeqRPM));
			  resOut.add(perRes);
			  perRes=null;		 
		}// for Rocket			
				 
		FileOperation.saveMatrixList(resOut, outFile);
		resOut=null;
  }
  
  List<SeqInfo> getQuerySeq(String querySeqFile){
		 
	 List<SeqInfo> querySeqList = new ArrayList<SeqInfo>();
	 SeqInfo seqInfo;
	 List<ArrayList<String>> queryStrList=FileOperation.getMatrixFromFile(querySeqFile);
	 for(int i=0;i<queryStrList.size();i++){
		 seqInfo= new SeqInfo();
		 seqInfo.seqName=queryStrList.get(i).get(0);
		 seqInfo.seq=queryStrList.get(i).get(1).toUpperCase();		
		 seqInfo.seqNumEncode=getSeqNumEncode(seqInfo.seq,maxEncodeBaseLen);
		 querySeqList.add(seqInfo);
		 seqInfo=null;
	 }
	 queryStrList=null;
	
	 
	 return querySeqList;
  }
  
  List<SeqRocket> buildQueryRockets(String querySeqFile){	 

		List<SeqInfo> querySeqList=getQuerySeq(querySeqFile);
		System.out.println("Total query objects: "+querySeqList.size());
		SeqRocket querySeqRocket;		
		SeqCompoRecognizer query;
		String querySeq="";		
		String rocketName="";		
		List<SeqRocket> rockets =new ArrayList<SeqRocket>();
		
		for(int i=0;i<querySeqList.size();i++){			
		    
		    if(querySeqList.get(i).seqName!=null)
			  rocketName=i+"_"+querySeqList.get(i).seqName;	
		    else rocketName="query_"+i;
			
		    querySeqRocket=new SeqRocket();
		    querySeqRocket.rocketName=rocketName;
		    querySeqRocket.isActive=true;
		
		    querySeqRocket.seqRecognizers=new ArrayList<SeqCompoRecognizer> ();
		    querySeqRocket.seqCompoFeatures=new SeqCompoFeatures();
		    querySeqRocket.seqCompoFeatures.compoNames=new ArrayList<String> ();	
		    querySeqRocket.seqCompoFeatures.compoColors=new ArrayList<String> ();
		    querySeqRocket.exactSeqCount=0;
		    querySeqRocket.nonExactSeqCount=0;
		    querySeqRocket.recognizedSeqCount=0;
		    
		    querySeq=querySeqList.get(i).seq;	
			if(querySeq!=null && querySeq.length()>0){			 		
			  querySeqRocket.seqCompoFeatures.compoNames.add(querySeqName);		
			  querySeqRocket.seqCompoFeatures.compoColors.add(querySeqColor);
			  
			  query=new SeqCompoRecognizer();	
			  query.index=querySeqRocket.seqCompoFeatures.compoNames.indexOf(querySeqName);				 
			  query.seq=querySeq;
			  query.rawSeq=querySeq;
			  query.seqName=rocketName+".QuerySeq";	
			  query.seqNumEncode=querySeqList.get(i).seqNumEncode;
		
			  querySeqRocket.seqRecognizers.add(query); 
			  query=null;
			}			
			
			rockets.add(querySeqRocket);
			querySeqRocket=null;		
		}
		
		configRecognizer(rockets);		
		return rockets;
  }
  

  void setNonExactAlignInfo(List<SeqRocket> querySeqRockets, String querySeqType,
		  List<SeqCompoAlignInfo> nonExactSeqList, int localExactEncodeLen){
		
	    if(querySeqRockets==null || nonExactSeqList==null) return;
	    
	    System.out.println(">>>Total reads of all NonExact matched query seq:"+nonExactSeqList.size());	
		//set seq encode
	    setSeqNumEncode(nonExactSeqList,localExactEncodeLen);
	    setSeqNumRevEncode(nonExactSeqList,localExactEncodeLen);
		Collections.sort(nonExactSeqList, new SeqCompoAlignInfo.CompSeqEncode(false));
		int compoIdx=querySeqRockets.get(0).seqCompoFeatures.compoNames.indexOf(querySeqName);
		String queryNonExactSeqFile;
		List<SeqCompoAlignInfo>queryNonExactSeqList;
		
		for(SeqRocket seqRocket:querySeqRockets){		
		   SeqCompoRecognizer querySeq=seqRocket.seqRecognizers.get(
			   seqRocket.seqCompoFeatures.compoNames.indexOf(querySeqType)
		   );
		   querySeq.seqNumEncode=getSeqNumEncode(querySeq.seq,localExactEncodeLen);
		   querySeq.seqNumRevEncode=getSeqNumRevEncode(querySeq.seq,localExactEncodeLen);
		   queryNonExactSeqList=findQuerySeq(nonExactSeqList,querySeq);		  
		   queryNonExactSeqList=getNoRecognizedSeq(queryNonExactSeqList,compoIdx);
	 	   querySeq=null;
		   System.out.println(seqRocket.rocketName+" NonExact reads:"+queryNonExactSeqList.size());
		   queryNonExactSeqFile=tmpDir+File.separator+"QueryNonExactMatchSeq_forward.fna";        
		   createBLASTTarSeq(queryNonExactSeqList,queryNonExactSeqFile);			
		   setNonExactBLASTInfo(seqRocket,queryNonExactSeqList,queryNonExactSeqFile,querySeqType);
		   queryNonExactSeqList=null;
		   FileOperation.delFile(queryNonExactSeqFile); 	
		  
		   System.out.println(seqRocket.rocketName+">> Got recognized seq:"
			      +" Exact:"+seqRocket.exactSeqCount  
			      +" NoExact:"+seqRocket.nonExactSeqCount);

		}// for seqRocket
  }
  
  void setNonExactBLASTInfo(SeqRocket seqRocket,List<SeqCompoAlignInfo>queryNonExactSeqList, 
		String queryNoExactSeqFile, String querySeqType){			
	      
	      if(!seqRocket.seqCompoFeatures.compoNames.contains(querySeqType)){
			  System.out.println("Warning: You did not set query sequence for '"
				       +seqRocket.rocketName+"', so no sequence extraction!!!");
			  
			  seqRocket.isDone=false;
			  return;
	      }
		  String tarBlastCMD="";
		  String blastCMD="";	
		  String querySeqFile="";		
		  String timeStamp=new SimpleDateFormat("yyyyMMdd_HHmmss").format(Calendar.getInstance().getTime());
		  String blastOutFile=tmpDir+File.separator+"Count.BlastOut."+timeStamp+".txt";
		  //List<Integer> blastWordSizeList;		 	
		  int blastWordSize=7;
		  String blastTask="blastn-short";	
		  int seqNum=queryNonExactSeqList.size();		
		  tarBlastCMD="-evalue 10000 -max_target_seqs "+seqNum;	
		  boolean doBlast=true;	
		  if(seqNum==0) doBlast=false;		  
		  
		  String seqType="";
		  SeqCompoRecognizer recognizer;	
		  List<SeqCompoAlignInfo> tarBlastSeqObj;	
		  
		  if(seqRocket.isActive){	
			for(int k=0;k<seqRocket.seqRecognizers.size();k++){			   
			   recognizer=seqRocket.seqRecognizers.get(k);
			   if(recognizer!=null){
				 seqType=seqRocket.seqCompoFeatures.compoNames.get(recognizer.index);
				 if(seqType.equals(querySeqType)){  	
				   tarBlastSeqObj=new ArrayList<SeqCompoAlignInfo>(); 			  
				   if(doBlast){
					  SeqOperation.createFASTASeq(recognizer.seq,recognizer.seqName,recognizer.seqFASTAFile); 
					  querySeqFile=recognizer.seqFASTAFile;
					  blastWordSize=recognizer.blastWordSize;
					  blastTask=recognizer.blastTask;			
					  blastCMD="blastn -task "+blastTask+" "
					          +tarBlastCMD+" -word_size="+blastWordSize
					          +" -query "+querySeqFile
					          +" -subject "+queryNoExactSeqFile
					          +" -out "+blastOutFile+" -outfmt 6";				
					  SeqOperation.runBLAST(blastCMD);			  
					  tarBlastSeqObj=getLeftSideBLASTSeq(queryNonExactSeqList,recognizer,
								  seqType,seqRocket,blastOutFile);		  
					  FileOperation.delFile(blastOutFile);
					  FileOperation.delFile(querySeqFile);
				   } 		              
				   seqRocket.nonExactSeqCount+=tarBlastSeqObj.size();
				   tarBlastSeqObj=null; //release memory space						 					 
				 }
				  
				 seqRocket.seqRecognizers.get(k).done=true;		
				 recognizer=null;	 
				 seqRocket.isDone=true;
			   }
			} // for seqRecognizer

		  }// is rocket active?         
    }
 
}
  
  
