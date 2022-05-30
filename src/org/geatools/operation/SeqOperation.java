package org.geatools.operation;
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
import java.io.*;
import java.nio.charset.StandardCharsets;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Calendar;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;

import org.geatools.data.structure.ChrSite;
import org.geatools.data.structure.SeqAlignSite;
import org.geatools.data.structure.SeqEncode;
import org.geatools.data.structure.SeqInfo;

public class SeqOperation{
 
   public static final String SEQTYPE_SINGLEEND="SingleEnd";
   public static final String SEQTYPE_PAIREND="PairEnd";
   public static int nameColIdx=0;
   public static int nameStartRowIdx=0;
   public static int splitStep=1000000;
 
   private static final String COMPLEMENT_TABLE 
	  // 0123456789ABCDEF0123456789ABCDEF
	  = "                                " // 0-31
	  + "             -                  " // 32-63
	  + " TVGH  CD  M KN   YSAABWXR      " // 64-95
	  + " tvgh  cd  m kn   ysaabwxr      "; // 96-127
	  //  ABCDEFGHIJKLMNOPQRSTUVWXYZ

   private static final byte[] COMPLEMENT_TABLE_BYTES 
	    = COMPLEMENT_TABLE.getBytes( StandardCharsets.US_ASCII );
 
   public SeqOperation(){

   }
	  
   public static byte[] reverseComplement( byte[] sequence ) {
      int length = sequence.length;
      byte[] result = new byte[ length ];

      for ( int i = 0; i < length; ++i ) {
        result[ (length - i) - 1] = COMPLEMENT_TABLE_BYTES[ sequence[i] ];
      }

      return result;
   }

   public static String reverseComplement( String sequence ) {
      byte[] complementBytes = reverseComplement(sequence.getBytes(StandardCharsets.US_ASCII));
      return new String(complementBytes, StandardCharsets.US_ASCII);
   } 
 
   public static byte[] complement( byte[] sequence ) {
      int length = sequence.length;
      byte[] result = new byte[ length ];

      for ( int i = 0; i < length; ++i ) {
        result[i] = COMPLEMENT_TABLE_BYTES[ sequence[i] ];
      }

      return result;
   }
	
   public static String complement( String sequence ) {
      byte[] complementBytes = complement(sequence.getBytes(StandardCharsets.US_ASCII));
      return new String(complementBytes, StandardCharsets.US_ASCII);
   }
 
   public static HashMap<Character, Integer> getBase2NumMap(){
	  
	  final HashMap<Character, Integer> base2numMap =new HashMap<Character, Integer>(); 
	  base2numMap.put('A', 1);
	  base2numMap.put('C', 2);
	  base2numMap.put('G', 3);
	  base2numMap.put('T', 4);
	  base2numMap.put('N', 0);
	  base2numMap.put('a', 1);
	  base2numMap.put('c', 2);
	  base2numMap.put('g', 3);
	  base2numMap.put('t', 4);
	  base2numMap.put('n', 0);
	  
	  return base2numMap;
   }
 
   public static long getSeqNumEncode(String seq, int encodeBaseLen){
	  if(encodeBaseLen<1) encodeBaseLen=1;
	  Integer val;
	  char ch;
	  long encodeVal=0;
	  String encodeValStr="";
	  int step=seq.length()/encodeBaseLen;
	  if(step==0) step=1;
	  int i=1;
	  for(int b = 0; b < seq.length(); b=b+step){
	    if(i>encodeBaseLen) break;
		i++;
		ch = seq.charAt(b);			   
		val =getBase2NumMap().get(ch);
		if(val!=null) encodeValStr=encodeValStr+val;		 
	  }
	  encodeVal=Long.parseLong(encodeValStr);
	 
	  return encodeVal;
   }
 
   public static long getSeqNumRevEncode(String seq,int encodeBaseLen){
	  if(encodeBaseLen<1) encodeBaseLen=1;
	  Integer val;
	  char ch;
	  long encodeVal=0;
	  String encodeValStr="";
	  int seqLen=seq.length();
	  int step=seq.length()/encodeBaseLen;
	  if(step==0) step=1;
	  int i=1;
	  for (int b = seqLen-1; b >=0; b=b-step) {	
		 if(i>encodeBaseLen) break;
		 i++;
		 ch = seq.charAt(b);			   
		 val =SeqOperation.getBase2NumMap().get(ch);
		 if(val!=null) encodeValStr=encodeValStr+val;			 
	  }
	  encodeVal=Long.parseLong(encodeValStr);
	 
	  return encodeVal;	
   }
 
   public static void setSeqNumEncode(List<SeqInfo> seqInfoList, float encodeBaseRate){
	  if(encodeBaseRate>1) encodeBaseRate=1;
	  if(encodeBaseRate<0) encodeBaseRate=0;
	  int encodeBaseLen;
	  for(SeqInfo seqInfo:seqInfoList){
	    encodeBaseLen=(int) (seqInfo.seq.length()*encodeBaseRate+0.5);
	    setSeqNumEncode(seqInfo,encodeBaseLen);
	  }	
   }
 
   public static void setSeqNumEncode(List<SeqInfo> seqInfoList, int encodeBaseLen){
      for(SeqInfo seqInfo:seqInfoList){
    	 setSeqNumEncode(seqInfo, encodeBaseLen);
      }	
   }
 
   public static void setSeqNumEncode(SeqInfo seqInfo, int encodeBaseLen){
	  if(encodeBaseLen<1) encodeBaseLen=1;
	  Integer val;
	  char ch;
	  long encodeVal=0;
	  String encodeValStr="";
	  int step=seqInfo.seq.length()/encodeBaseLen;
	  if(step==0) step=1;
	  int i=1;
	  for(int b = 0; b < seqInfo.seq.length(); b=b+step){
		if(i>encodeBaseLen) break;
		i++;
		ch = seqInfo.seq.charAt(b);			   
		val =getBase2NumMap().get(ch);
		if(val!=null) encodeValStr=encodeValStr+val;		 
	  }
	  encodeVal=Long.parseLong(encodeValStr);
	  seqInfo.seqNumEncode=encodeVal;	
   }
 
 /*
   public static void setSeqNumEncode(List<SeqInfo> seqInfoList, int [] baseIndices){
      for(SeqInfo seqInfo:seqInfoList){
    	 setSeqEncodeVal(seqInfo, baseIndices);
      }	
   }
 
   public static void setSeqNumEncode(SeqInfo seqInfo, int [] baseIndices){
	  Integer val;
	  char ch;
	  long encodeVal=0;
	  String encodeValStr="";
	  int b;
	  for (int i = 0; i < baseIndices.length; i++) {
		 b=baseIndices[i];
		 if(b>(seqInfo.seq.length()-1)) b=seqInfo.seq.length()-1;
		 ch = seqInfo.seq.charAt(b);			   
		 val =getBase2NumMap().get(ch);
		 if(val!=null)	encodeValStr=encodeValStr+val;
		 
	  }
	  encodeVal=Long.parseLong(encodeValStr);
	  seqInfo.seqNumEncode=encodeVal;	
   }
 */
 
   public static void setSeqNumRevEncode(List<SeqInfo> seqInfoList, float encodeBaseRate){
	  if(encodeBaseRate>1) encodeBaseRate=1;
	  if(encodeBaseRate<0) encodeBaseRate=0;
	  int encodeBaseLen;
	  for(SeqInfo seqInfo:seqInfoList){
	   	encodeBaseLen=(int) (seqInfo.seq.length()*encodeBaseRate+0.5);
	   	setSeqNumRevEncode(seqInfo,encodeBaseLen);
	  }	
   }
 
   public static void setSeqNumRevEncode(List<SeqInfo> seqInfoList, int encodeBaseLen){
      for(SeqInfo seqInfo:seqInfoList){
    	  setSeqNumRevEncode(seqInfo,encodeBaseLen);
      }	
   }	    

   public static void setSeqNumRevEncode(SeqInfo seqInfo,int encodeBaseLen){
	  if(encodeBaseLen<1) encodeBaseLen=1;
	  Integer val;
	  char ch;
	  long encodeVal=0;
	  String encodeValStr="";
	  int seqLen=seqInfo.seq.length();
	  int step=seqInfo.seq.length()/encodeBaseLen;
	  if(step==0) step=1;
	  int i=1;
	  for (int b = seqLen-1; b >=0; b=b-step) {	
		 if(i>encodeBaseLen) break;
		 i++;
		 ch = seqInfo.seq.charAt(b);			   
		 val =SeqOperation.getBase2NumMap().get(ch);
		 if(val!=null) encodeValStr=encodeValStr+val;			 
	  }
	  encodeVal=Long.parseLong(encodeValStr);
	  seqInfo.seqNumRevEncode=encodeVal;		
   }
  
   public static List<SeqEncode> getSeqNumEncode(String inSeqFile, int splitStep, int encodeBaseLen){
	     
	  List<SeqEncode> seqEnocdeObjList=new ArrayList<SeqEncode>();	
	     
	  if(isFASTASeq(inSeqFile)){			
	   	seqEnocdeObjList=getNumEncodeFromFASTA(inSeqFile, splitStep, encodeBaseLen);		
	  }else if(isFASTQSeq(inSeqFile)){			 
		seqEnocdeObjList=getNumEncodeFromFASTQ(inSeqFile, splitStep, encodeBaseLen);	
	  }
	     
	  return seqEnocdeObjList;
   }

   public static List<SeqEncode> getNumEncodeFromFASTA(String inSeqFile,int splitStep,int encodeBaseLen){
	    	   
	    List<SeqEncode> seqEnocdeObjList=new ArrayList<SeqEncode>();	   
	    List<SeqInfo> seqObjList;	
		if(splitStep<0) splitStep=500000;

		int total=SeqOperation.getFASTASeqNum(inSeqFile);
		System.out.println("Total reads to be encoded: "+total); 
		int stopPos=0;		 
		for(int f=0;f<total;f=f+splitStep){	      
	        stopPos=Math.min(total, (f+splitStep));
	        seqObjList=getFASTASeqObj(inSeqFile,f,stopPos);	
	        seqEnocdeObjList.addAll(getSeqNumEncode(seqObjList,encodeBaseLen));
	        seqObjList=null;
		    System.out.println((f+1)+"-"+stopPos+" encoded!");	
		}
		  
		return  seqEnocdeObjList;
   }

   public static List<SeqEncode> getNumEncodeFromFASTQ(String inSeqFile,int splitStep,int encodeBaseLen){
	   
	    List<SeqEncode> seqEnocdeObjList=new ArrayList<SeqEncode>();	   
	    List<SeqInfo> seqObjList;	
		if(splitStep<0) splitStep=1000000;

		int total=SeqOperation.getFASTQSeqNum(inSeqFile);
		System.out.println("Total reads to be encoded: "+total); 
		int stopPos=0;		 
		for(int f=0;f<total;f=f+splitStep){	      
	        stopPos=Math.min(total, (f+splitStep));
	        seqObjList=getFASTQSeqObj(inSeqFile,f,stopPos);	
	        seqEnocdeObjList.addAll(getSeqNumEncode(seqObjList,encodeBaseLen));
	        seqObjList=null;
		    System.out.println((f+1)+"-"+stopPos+" encoded!");	
		}
		  
		return  seqEnocdeObjList;
   }
 
   public static List<SeqEncode> getSeqNumEncode(List<SeqInfo> seqObjList, int encodeBaseLen){
	 
	  List<SeqEncode> seqEnocdeObjList=new ArrayList<SeqEncode>();
	 
	  for(SeqInfo seqObj:seqObjList){
		 SeqEncode perSeq=new SeqEncode();
		 perSeq.seqID=seqObj.seqID;
		 perSeq.numEncode=getSeqNumEncode(seqObj.seq,encodeBaseLen);
		 perSeq.numRevEncode=getSeqNumRevEncode(seqObj.seq,encodeBaseLen);	
		 seqEnocdeObjList.add(perSeq);
		 perSeq=null;
      }	
	 
	  return seqEnocdeObjList;
   }	
  
   public static long getMaxSeqNumEncode(List<SeqEncode> seqEncodeObjList){
	  long num=Long.MIN_VALUE;
	  
	  for(SeqEncode seq:seqEncodeObjList) {
		 if(seq.numEncode>num) num=seq.numEncode;
	  }
	  
	  return num;
   }
  
   public static long getMinSeqNumEncode(List<SeqEncode> seqEncodeObjList){
	  long num=Long.MAX_VALUE;
	  
	  for(SeqEncode seq:seqEncodeObjList) {
		 if(seq.numEncode<num) num=seq.numEncode;
	  }
	  
	  return num;
   }
   
   public static List<String> getEncodeInterval(List<SeqEncode> seqEncodeObjList, int intervalNum){	  
	  
	  List<String> outList=new ArrayList<String>();
	  boolean desc=false;
	  Collections.sort(seqEncodeObjList, new SeqEncode.CompSeqEncode(desc));
	  long num1 = 0;
	  long num2; 
	  int sizePerInterval=(int) Math.ceil((1.0d*seqEncodeObjList.size())/(1.0d*intervalNum));
	  int k=1;
	  int s=1;
	  for(int i=0;i<seqEncodeObjList.size();i++) {
		 if(k==1) num1=seqEncodeObjList.get(i).numEncode;		 
		 if(k>=sizePerInterval || i==seqEncodeObjList.size()-1){
			num2=seqEncodeObjList.get(i).numEncode;
			if(i==seqEncodeObjList.size()-1) {
			   outList.add(num1+";"+num2+";encode_interval"+s+"."+k);
			}else if(num2<seqEncodeObjList.get(i+1).numEncode){
			   outList.add(num1+";"+num2+";encode_interval"+s+"."+k);
			   k=1;
			   s++;
			}else {
			   k++;
			}			
		 }else {
			k++; 
		 }		
	  }
	  
	  return outList;
   }
    
   public static List<String> getEncodeInterval0(long minNumEncode,long maxNumEncode,long numEncodeInterval){
	  long num1;
	  long num2;
	  List<String> outList=new ArrayList<String>();
	  int k=1;
	  for(long i=minNumEncode;i<=maxNumEncode;i=i+numEncodeInterval) {
		  num1=i;
		  num2=i+numEncodeInterval;
		  outList.add(num1+";"+num2+";encode_interval"+k);
		  k++;
	  }
	  return outList;
   }
  
   public static List<String> dispatchSeqByNumEncode(String inSeqFile,int splitStep, 
		  List<SeqEncode> seqEncodeObjList,int splitNum, String outDir){
	     
	  List<String> dividedFiles=new ArrayList<String>();	 	     
	    
	  if(isFASTASeq(inSeqFile)){			
		  dividedFiles=dispatchFASTAByNumEncode(inSeqFile,splitStep,seqEncodeObjList,splitNum,outDir);		
	  }else if(isFASTQSeq(inSeqFile)){			 
		  dividedFiles=dispatchFASTQByNumEncode(inSeqFile,splitStep,seqEncodeObjList,splitNum,outDir);
	  }
	     
	  return dividedFiles;
   }

   public static List<String> dispatchFASTAByNumEncode(String inSeqFile,int seqSplitStep, 
		  List<SeqEncode> seqEncodeObjList,int intervalNum, String outDir){
	    		
		List<String> encodeIntervalList=getEncodeInterval(seqEncodeObjList,intervalNum);
		
		boolean desc=false;
		Collections.sort(seqEncodeObjList, new SeqEncode.CompSeqID(desc));
		
		List<String> dividedFiles=new ArrayList<String>();
	    List<SeqInfo> seqObjList;	
		if(outDir==null) outDir=inSeqFile.substring(0,inSeqFile.lastIndexOf("/"))+"/"
		                   +inSeqFile.substring(inSeqFile.lastIndexOf("/")+1,inSeqFile.lastIndexOf("."))
		                   +"_EncodeDivided";
		FileOperation.newFolder(outDir);
		int total=SeqOperation.getFASTASeqNum(inSeqFile);
		System.out.println("Total reads to be dispatched: "+total);  		
		
		try{
		  String dividedFileTag=inSeqFile.substring(inSeqFile.lastIndexOf("/")+1,inSeqFile.lastIndexOf("."));
		  List<Long> numEncode1List=new ArrayList<Long>();
		  List<Long> numEncode2List=new ArrayList<Long>();
		  List<BufferedWriter> writerList=new ArrayList<BufferedWriter>();
		  List<BufferedWriter> writerEncodeList=new ArrayList<BufferedWriter>();
		  BufferedWriter writer=null;
		  String[] items;

		  for(String str:encodeIntervalList) {
			items=str.split(";");
			if(items.length>2) {
			  numEncode1List.add(Long.parseLong(items[0]));
			  numEncode2List.add(Long.parseLong(items[1]));
			  String outFile=outDir+"/"+dividedFileTag+"."+items[2]+".fna";			 
			  writer=new BufferedWriter(new FileWriter(outFile));
			  writerList.add(writer);
			  writer=null;
			  dividedFiles.add(outFile);
			  outFile=null;
			  
			  outFile=outDir+"/"+dividedFileTag+"."+items[2]+".encode";			 
			  writer=new BufferedWriter(new FileWriter(outFile));
			  writerEncodeList.add(writer);
			  writer=null;		
			  outFile=null;			  
			}
		  }

		  int stopPos;	
		  int seqID;
		  String seqIdentifier;
		  String seqName;
		  String seq;	
		  long numEncode;
		  long numRevEncode;	 
		  for(int f=0;f<total;f=f+seqSplitStep){	         
	         stopPos=Math.min(total, (f+seqSplitStep));
	         seqObjList=getFASTASeqObj(inSeqFile,f,stopPos);	
	         for(int i=0;i<seqObjList.size();i++){	
	        	seqID=seqObjList.get(i).seqID;
	        	seqIdentifier=seqObjList.get(i).seqIdentifier;
	        	seqName=seqObjList.get(i).seqName;
		 		seq=seqObjList.get(i).seq;
				numEncode=seqEncodeObjList.get(seqID).numEncode;
				numRevEncode=seqEncodeObjList.get(seqID).numRevEncode;
		 		writer=null;
		 		for(int e=0;e<numEncode1List.size();e++) {
		 		   if(numEncode>=numEncode1List.get(e) && numEncode<=numEncode2List.get(e)) {
				 	  writerList.get(e).write(">"+seqIdentifier);
				 	  writerList.get(e).newLine();				 
				 	  writerList.get(e).write(seq);
				 	  writerList.get(e).newLine();
				 	  writerList.get(e).flush();	
				 	  
				 	  writerEncodeList.get(e).write(seqName+"\t"+numEncode+"\t"+numRevEncode);
				 	  //writerEncodeList.get(e).write(seqName+"\t"+numEncode);
				 	  writerEncodeList.get(e).newLine();				
				 	  writerEncodeList.get(e).flush();	
				 	  
		 			  break;
		 		  }
		 		}			
		     }
	         seqObjList=null;
			 seqIdentifier=null;
			 seqName=null;
		     seq=null;
		     System.out.println((f+1)+"-"+stopPos+" dispatched!");	
		  }
		  for(int w=0;w<writerList.size();w++) {
			  writerList.get(w).close();
			  writerEncodeList.get(w).close();
		  }
		  writerList=null;
		  writerEncodeList=null;
		}catch(IOException e){
	        System.out.println(e);
	    }
		
		encodeIntervalList=null;
		  
		return dividedFiles;
   }

   public static List<String> dispatchFASTQByNumEncode(String inSeqFile,int seqSplitStep,
		  List<SeqEncode> seqEncodeObjList, int intervalNum, String outDir){  
		
		List<String> encodeIntervalList=getEncodeInterval(seqEncodeObjList,intervalNum);
		
		boolean desc=false;
		Collections.sort(seqEncodeObjList, new SeqEncode.CompSeqID(desc));
		
		List<String> dividedFiles=new ArrayList<String>();
	    List<SeqInfo> seqObjList;	
		if(outDir==null) outDir=inSeqFile.substring(0,inSeqFile.lastIndexOf("/"))+"/"
		                   +inSeqFile.substring(inSeqFile.lastIndexOf("/")+1,inSeqFile.lastIndexOf("."))
		                   +"_EncodeDivided";
		FileOperation.newFolder(outDir);
		int total=SeqOperation.getFASTQSeqNum(inSeqFile);
		System.out.println("Total reads to be dispatched: "+total);  		
		
		try{
		  String dividedFileTag=inSeqFile.substring(inSeqFile.lastIndexOf("/")+1,inSeqFile.lastIndexOf("."));
		  List<Long> numEncode1List=new ArrayList<Long>();
		  List<Long> numEncode2List=new ArrayList<Long>();
		  List<BufferedWriter> writerList=new ArrayList<BufferedWriter>();
		  List<BufferedWriter> writerEncodeList=new ArrayList<BufferedWriter>();
		  BufferedWriter writer=null;
		  String[] items;

		  for(String str:encodeIntervalList) {
			items=str.split(";");
			if(items.length>2) {
			  numEncode1List.add(Long.parseLong(items[0]));
			  numEncode2List.add(Long.parseLong(items[1]));
			  String outFile=outDir+"/"+dividedFileTag+"."+items[2]+".fastq";			 
			  writer=new BufferedWriter(new FileWriter(outFile));
			  writerList.add(writer);
			  writer=null;
			  dividedFiles.add(outFile);
			  outFile=null;
			  
			  outFile=outDir+"/"+dividedFileTag+"."+items[2]+".encode";			 
			  writer=new BufferedWriter(new FileWriter(outFile));
			  writerEncodeList.add(writer);
			  writer=null;		
			  outFile=null;		
			}
		  }

		  int stopPos;	
		  int seqID;
		  String seqIdentifier;
		  String seqName;
		  String seq;
		  String seqQuality;
		  long numEncode;
		  long numRevEncode;		  
		  for(int f=0;f<total;f=f+seqSplitStep){   
		     stopPos=Math.min(total, (f+seqSplitStep));
		     seqObjList=getFASTQSeqObj(inSeqFile,f,stopPos);	
		     for(int i=0;i<seqObjList.size();i++){	
		    	seqID=seqObjList.get(i).seqID;
		    	seqIdentifier=seqObjList.get(i).seqIdentifier;
		    	seqName=seqObjList.get(i).seqName;
				seq=seqObjList.get(i).seq;
				seqQuality=seqObjList.get(i).seqQualityEncode;
				numEncode=seqEncodeObjList.get(seqID).numEncode;
				numRevEncode=seqEncodeObjList.get(seqID).numRevEncode;

				writer=null;
				for(int e=0;e<numEncode1List.size();e++) {
				  if(numEncode>=numEncode1List.get(e) && numEncode<=numEncode2List.get(e)) {
					 writerList.get(e).write("@"+seqIdentifier);
					 writerList.get(e).newLine();		
					 writerList.get(e).write(seq);
					 writerList.get(e).newLine();
					 writerList.get(e).write("+");
					 writerList.get(e).newLine();
					 writerList.get(e).write(seqQuality);
					 writerList.get(e).newLine();
					 writerList.get(e).flush();	
					
				 	 writerEncodeList.get(e).write(seqName+"\t"+numEncode+"\t"+numRevEncode);
				 	 //writerEncodeList.get(e).write(seqName+"\t"+numEncode);
				 	 writerEncodeList.get(e).newLine();				
				 	 writerEncodeList.get(e).flush();	
				 	  
		 			 break;
		 		  }
		 		}
			
		     }
	         seqObjList=null;		     
			 seqIdentifier=null;
			 seqName=null;
		     seq=null;
		     seqQuality=null;
   
		     System.out.println((f+1)+"-"+stopPos+" dispatched!");	
		  }
		  for(int w=0;w<writerList.size();w++) {
			 writerList.get(w).close();
			 writerEncodeList.get(w).close();
		  }
		  writerList=null;
		  writerEncodeList=null;
		}catch(IOException e){
	        System.out.println(e);
	    }
			
		encodeIntervalList=null;
	    
		return dividedFiles;
   }
   
   
   public static List<ArrayList<String>> dispatchSeqByNumEncode(String inSeqFile,String inSeqFile2,
		      int splitStep, List<SeqEncode> seqEncodeObjList,int splitNum,String outDir){
		     
	    List<ArrayList<String>> dividedFiles=new ArrayList<ArrayList<String>>();	 	     
		    
		if(isFASTASeq(inSeqFile)){			
		    dividedFiles=dispatchFASTAByNumEncode(inSeqFile,inSeqFile2,splitStep,seqEncodeObjList,
					  splitNum,outDir);		
		}else if(isFASTQSeq(inSeqFile)){			 
		    dividedFiles=dispatchFASTQByNumEncode(inSeqFile,inSeqFile2,splitStep,seqEncodeObjList,
					  splitNum,outDir);
		}
		     
		return dividedFiles;
   }

   public static List<ArrayList<String>> dispatchFASTAByNumEncode(String inSeqFile, String inSeqFile2,
		      int seqSplitStep, List<SeqEncode> seqEncodeObjList, int intervalNum, String outDir){
					
		List<String> encodeIntervalList=getEncodeInterval(seqEncodeObjList,intervalNum);
	    
	    boolean desc=false;
	    Collections.sort(seqEncodeObjList, new SeqEncode.CompSeqID(desc));
			
		ArrayList<String> dividedFiles=new ArrayList<String>();
		ArrayList<String> dividedFiles2=new ArrayList<String>();
		List<SeqInfo> seqObjList;	
		List<SeqInfo> seqObjList2;
		if(outDir==null) outDir=inSeqFile.substring(0,inSeqFile.lastIndexOf("/"))+"/"
			                   +inSeqFile.substring(inSeqFile.lastIndexOf("/")+1,inSeqFile.lastIndexOf("."))
			                   +"_EncodeDivided";
		FileOperation.newFolder(outDir);
		int total=SeqOperation.getFASTASeqNum(inSeqFile);
		System.out.println("Total reads to be dispatched: "+total);  		
			
		try{
			String dividedFileTag=inSeqFile.substring(inSeqFile.lastIndexOf("/")+1,inSeqFile.lastIndexOf("."));
			String dividedFileTag2=inSeqFile2.substring(inSeqFile2.lastIndexOf("/")+1,inSeqFile2.lastIndexOf("."));
			List<Long> numEncode1List=new ArrayList<Long>();
			List<Long> numEncode2List=new ArrayList<Long>();
			List<BufferedWriter> writerList=new ArrayList<BufferedWriter>();
			List<BufferedWriter> writerEncodeList=new ArrayList<BufferedWriter>();
			List<BufferedWriter> writerList2=new ArrayList<BufferedWriter>();		
			BufferedWriter writer=null;
			String[] items;

			for(String str:encodeIntervalList) {
			   items=str.split(";");
			   if(items.length>2) {
				  numEncode1List.add(Long.parseLong(items[0]));
				  numEncode2List.add(Long.parseLong(items[1]));
				  //forward
				  String outFile=outDir+"/"+dividedFileTag+"."+items[2]+".fna";			 
				  writer=new BufferedWriter(new FileWriter(outFile));
				  writerList.add(writer);
				  writer=null;
				  dividedFiles.add(outFile);
				  outFile=null;
				  
				  outFile=outDir+"/"+dividedFileTag+"."+items[2]+".encode";			 
				  writer=new BufferedWriter(new FileWriter(outFile));
				  writerEncodeList.add(writer);
				  writer=null;		
				  outFile=null;	
				  
				  //reverse
				  outFile=outDir+"/"+dividedFileTag2+"."+items[2]+".fna";			 
				  writer=new BufferedWriter(new FileWriter(outFile));
				  writerList2.add(writer);
				  writer=null;
				  dividedFiles2.add(outFile);
				  outFile=null;	 
			   }
			}
			  
			int seqID;
			String seqIdentifier;
			String seqName;
			String seq;
			long numEncode;
			long numRevEncode;
			  
			//int seqID2;
			String seqIdentifier2;
			String seq2;
			  
			int stopPos;
			for(int f=0;f<total;f=f+seqSplitStep){	         
		        stopPos=Math.min(total, (f+seqSplitStep));
		        seqObjList=getFASTASeqObj(inSeqFile,f,stopPos);	
		        seqObjList2=getFASTASeqObj(inSeqFile2,f,stopPos);
		        for(int i=0;i<seqObjList.size();i++){	
		        	seqID=seqObjList.get(i).seqID;
		        	seqIdentifier=seqObjList.get(i).seqIdentifier;
		        	seqName=seqObjList.get(i).seqName;
			 		seq=seqObjList.get(i).seq;
					numEncode=seqEncodeObjList.get(seqID).numEncode;
					numRevEncode=seqEncodeObjList.get(seqID).numRevEncode;
					
		        	//seqID2=seqObjList2.get(i).seqID;
		        	seqIdentifier2=seqObjList2.get(i).seqIdentifier;
			 		seq2=seqObjList2.get(i).seq;
					
			 		for(int e=0;e<numEncode1List.size();e++) {
			 		   if(numEncode>=numEncode1List.get(e) && numEncode<=numEncode2List.get(e)) {
					 	  //forward
			 			  writerList.get(e).write(">"+seqIdentifier);
					 	  writerList.get(e).newLine();					  
					 	  writerList.get(e).write(seq);
					 	  writerList.get(e).newLine();
					 	  writerList.get(e).flush();	
					 	  
					 	  writerEncodeList.get(e).write(seqName+"\t"+numEncode+"\t"+numRevEncode);
					 	  //writerEncodeList.get(e).write(seqName+"\t"+numEncode);
					 	  writerEncodeList.get(e).newLine();				
					 	  writerEncodeList.get(e).flush();	
					 	  
					 	  //reverse
			 			  writerList2.get(e).write(">"+seqIdentifier2);
					 	  writerList2.get(e).newLine();					 	  
					 	  writerList2.get(e).write(seq2);
					 	  writerList2.get(e).newLine();
					 	  writerList2.get(e).flush();	
					 	  
			 			  break;
			 		  }
			 		}			
			    }
		        seqObjList=null;
				seqIdentifier=null;
				seqName=null;
			    seq=null;
		        seqObjList2=null;
				seqIdentifier2=null;
			    seq2=null;
			    System.out.println((f+1)+"-"+stopPos+" dispatched!");	
			}
			for(int w=0;w<writerList.size();w++) {
				 writerList.get(w).close();
				 writerEncodeList.get(w).close();
				 writerList2.get(w).close();		
			}
			writerList=null;
			writerEncodeList=null;
			writerList2=null;	
		}catch(IOException e){
		    System.out.println(e);
		}
			
		encodeIntervalList=null;
			
		List<ArrayList<String>> dividedFiles12=new ArrayList<ArrayList<String>>();
		dividedFiles12.add(dividedFiles);
		dividedFiles12.add(dividedFiles2);
			
		return dividedFiles12;
   }

   public static List<ArrayList<String>> dispatchFASTQByNumEncode(String inSeqFile, String inSeqFile2,
		     int seqSplitStep, List<SeqEncode> seqEncodeObjList,int intervalNum,String outDir){	    
		  
			
		List<String> encodeIntervalList=getEncodeInterval(seqEncodeObjList,intervalNum);
		
		boolean desc=false;
		Collections.sort(seqEncodeObjList, new SeqEncode.CompSeqID(desc));
			
		ArrayList<String> dividedFiles=new ArrayList<String>();
		ArrayList<String> dividedFiles2=new ArrayList<String>();
		List<SeqInfo> seqObjList;	
		List<SeqInfo> seqObjList2;
		  
		if(outDir==null) outDir=inSeqFile.substring(0,inSeqFile.lastIndexOf("/"))+"/"
			                   +inSeqFile.substring(inSeqFile.lastIndexOf("/")+1,inSeqFile.lastIndexOf("."))
			                   +"_EncodeDivided";
		FileOperation.newFolder(outDir);
		int total=SeqOperation.getFASTQSeqNum(inSeqFile);
		System.out.println("Total reads to be dispatched: "+total);  		
			
		try{
			String dividedFileTag=inSeqFile.substring(inSeqFile.lastIndexOf("/")+1,inSeqFile.lastIndexOf("."));
			String dividedFileTag2=inSeqFile2.substring(inSeqFile2.lastIndexOf("/")+1,inSeqFile2.lastIndexOf("."));
			List<Long> numEncode1List=new ArrayList<Long>();
			List<Long> numEncode2List=new ArrayList<Long>();
			List<BufferedWriter> writerList=new ArrayList<BufferedWriter>();
			List<BufferedWriter> writerEncodeList=new ArrayList<BufferedWriter>();
			List<BufferedWriter> writerList2=new ArrayList<BufferedWriter>();		
			BufferedWriter writer=null;
			String[] items;

			for(String str:encodeIntervalList) {
			   items=str.split(";");
			   if(items.length>2) {
				  numEncode1List.add(Long.parseLong(items[0]));
				  numEncode2List.add(Long.parseLong(items[1]));
				  
				  //forward
				  String outFile=outDir+"/"+dividedFileTag+"."+items[2]+".fastq";			 
				  writer=new BufferedWriter(new FileWriter(outFile));
				  writerList.add(writer);
				  writer=null;
				  dividedFiles.add(outFile);
				  outFile=null;
				  
				  outFile=outDir+"/"+dividedFileTag+"."+items[2]+".encode";			 
				  writer=new BufferedWriter(new FileWriter(outFile));
				  writerEncodeList.add(writer);
				  writer=null;		
				  outFile=null;	
				  
				  //reverse
				  outFile=outDir+"/"+dividedFileTag2+"."+items[2]+".fastq";			 
				  writer=new BufferedWriter(new FileWriter(outFile));
				  writerList2.add(writer);
				  writer=null;
				  dividedFiles2.add(outFile);
				  outFile=null;
			   }
			}			  	
			  
			int seqID;
			String seqIdentifier;
			String seqName;
			String seq;
			String seqQuality;
			long numEncode;
			long numRevEncode;	
			  
			//int seqID2;
			String seqIdentifier2;
			String seq2;
			String seqQuality2;
			  
	        int stopPos;
			for(int f=0;f<total;f=f+seqSplitStep){   
			   stopPos=Math.min(total, (f+seqSplitStep));
			   seqObjList=getFASTQSeqObj(inSeqFile,f,stopPos);
			   seqObjList2=getFASTQSeqObj(inSeqFile2,f,stopPos);	
			   for(int i=0;i<seqObjList.size();i++){	
			      seqID=seqObjList.get(i).seqID;
			      seqIdentifier=seqObjList.get(i).seqIdentifier;
			      seqName=seqObjList.get(i).seqName;
				  seq=seqObjList.get(i).seq;
				  seqQuality=seqObjList.get(i).seqQualityEncode;
				  numEncode=seqEncodeObjList.get(seqID).numEncode;
				  numRevEncode=seqEncodeObjList.get(seqID).numRevEncode;
					
		          //seqID2=seqObjList2.get(i).seqID;
		          seqIdentifier2=seqObjList2.get(i).seqIdentifier;
			 	  seq2=seqObjList2.get(i).seq;
			 	  seqQuality2=seqObjList2.get(i).seqQualityEncode;

				  for(int e=0;e<numEncode1List.size();e++) {
					  if(numEncode>=numEncode1List.get(e) && numEncode<=numEncode2List.get(e)) {
						 //forward
						 writerList.get(e).write("@"+seqIdentifier);
						 writerList.get(e).newLine();		
						 writerList.get(e).write(seq);
						 writerList.get(e).newLine();
						 writerList.get(e).write("+");
						 writerList.get(e).newLine();
						 writerList.get(e).write(seqQuality);
						 writerList.get(e).newLine();
						 writerList.get(e).flush();	
						
					 	 writerEncodeList.get(e).write(seqName+"\t"+numEncode+"\t"+numRevEncode);
						 //writerEncodeList.get(e).write(seqName+"\t"+numEncode);
						 writerEncodeList.get(e).newLine();				
					 	 writerEncodeList.get(e).flush();	
					 	
						 //reverse
						 writerList2.get(e).write("@"+seqIdentifier2);
						 writerList2.get(e).newLine();		
						 writerList2.get(e).write(seq2);
						 writerList2.get(e).newLine();
						 writerList2.get(e).write("+");
						 writerList2.get(e).newLine();
						 writerList2.get(e).write(seqQuality2);
						 writerList2.get(e).newLine();
						 writerList2.get(e).flush();	
					 	  
			 			 break;
			 		  }
			 	  }				
			   }
		       seqObjList=null;		     
			   seqIdentifier=null;
			   seqName=null;
			   seq=null;
			   seqQuality=null;
			     
		       seqObjList2=null;		     
			   seqIdentifier2=null;
			   seq2=null;
			   seqQuality2=null;
	   
			   System.out.println((f+1)+"-"+stopPos+" dispatched!");	
			}
			for(int w=0;w<writerList.size();w++) {
			   writerList.get(w).close();
			   writerEncodeList.get(w).close();
			   writerList2.get(w).close();
			}
			writerList=null;
			writerEncodeList=null;
			writerList2=null;		
		}catch(IOException e){
		    System.out.println(e);
		}
				
		encodeIntervalList=null;
		    
		List<ArrayList<String>> dividedFiles12=new ArrayList<ArrayList<String>>();
		dividedFiles12.add(dividedFiles);
		dividedFiles12.add(dividedFiles2);
				
		return dividedFiles12;
   }

 /*
   public static List<SeqEncode> findSeqEncodeObj(List<SeqEncode>sortedSeqEncodeList, SeqEncode querySeq){
		 
	     List<SeqEncode> querySeqList= new ArrayList<SeqEncode>();
	     int seqIdx = Collections.binarySearch(
	    		 sortedSeqEncodeList,
				 querySeq, 
				 new SeqEncode.CompSeqEncode(false)
		 );

		 if(seqIdx>=0){		  
			 //go backward......
			 for(int j=seqIdx;j>=0;j--){			   
			   if(querySeq.numEncode!=sortedSeqEncodeList.get(j).numEncode) break;
			  
			   if(querySeq.numEncode==sortedSeqEncodeList.get(j).numEncode)
			     querySeqList.add(sortedSeqEncodeList.get(j));
			 }
			 //go forward......
		     for(int j=seqIdx+1;j<sortedSeqEncodeList.size();j++){		     
		       if(querySeq.numEncode!=sortedSeqEncodeList.get(j).numEncode) break;
			   
			   if(querySeq.numEncode==sortedSeqEncodeList.get(j).numEncode)
		         querySeqList.add(sortedSeqEncodeList.get(j));
		     }
		 }
		 
		 return querySeqList;
   }
*/
   public static List<SeqInfo> findSeqObjOfSameEncode(List<SeqInfo> sortedSeqList,SeqInfo querySeq){
	 
       List<SeqInfo> querySeqList= new ArrayList<SeqInfo>();
       int seqIdx = Collections.binarySearch(
    		 sortedSeqList,
			 querySeq, 
			 new SeqInfo.CompSeqEncode(false)
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
 
   public static int runBLAST(String blastCMD){
	   try {
	     Runtime rt = Runtime.getRuntime();        
	     Process pr = rt.exec(blastCMD);
	 
	     BufferedReader input = new BufferedReader(new InputStreamReader(pr.getInputStream()));
	     String line=null;
	     while((line=input.readLine()) != null) {
	       System.out.println(line);
	     }
	 
	     int exitVal = pr.waitFor();
		 if(exitVal>0) System.out.println("Error: Exited with error code "+exitVal);
		 pr=null;
			
		 return exitVal;	 
	  }catch(Exception e) {
	     System.out.println(e.toString());
	     e.printStackTrace();
	        
	     return 1;
	  }	    
   }
 
   public static void createFASTASeq(String seq, String seqName, String outSeqFile){
		  
	  try{    
	    BufferedWriter writer=null;
	    writer=new BufferedWriter(new FileWriter(outSeqFile));
	    String nameLine;
		//int seqID=0;
		
	    nameLine=">"+seqName+" "+seq.length();
		writer.write(nameLine);
		writer.newLine();
		  //writer.flush(); 
		writer.write(seq);
		writer.newLine();
		writer.flush(); 
			
		writer.close();
			
	  }catch(IOException e){
	    System.out.println(e);
	  }  
	 
   }
 
   public static boolean isFASTASeq(List<String> seqFileList){
	   if(seqFileList==null || seqFileList.size()==0) return false;
	   boolean isOK=false;
	   for(String seqFile:seqFileList){
		 if(isFASTASeq(seqFile)){
			 isOK=true;
			 break;
		 }
	   }
	 
	   return isOK;
   }
 
   public static boolean isFASTASeq(String seqFile){
	   if(seqFile==null) return false;
	 
	   File f=new File(seqFile);
	   if(!f.exists()){
		  f=null; 
		  return false;
	   }
	   
	   boolean isOk=false;
	 
	   String inSeqFileFormat=FileOperation.getFileFormat(seqFile);
	
       if(inSeqFileFormat.equalsIgnoreCase("fasta") 
			|| inSeqFileFormat.equalsIgnoreCase("fna") 
			|| inSeqFileFormat.equalsIgnoreCase("fa")){
		
		 isOk=true;	
	   }
	 
	   return isOk;
   }
 
   public static boolean isFASTQSeq(List<String> seqFileList){
	   if(seqFileList==null || seqFileList.size()==0) return false;
	   boolean isOK=false;
	   for(String seqFile:seqFileList){
		 if(isFASTQSeq(seqFile)){
			 isOK=true;
			 break;
		 }
	   }
	 
	   return isOK;
   }
 
   public static boolean isFASTQSeq(String seqFile){
	   if(seqFile==null) return false;
	 
	   File f=new File(seqFile);
	   if(!f.exists()){
		 f=null; 
		 return false;
	   }
	 
	   boolean isOk=false;
	 
	   String inSeqFileFormat=FileOperation.getFileFormat(seqFile);
	
       if(inSeqFileFormat.equalsIgnoreCase("fastq") 
			|| inSeqFileFormat.equalsIgnoreCase("fnq") 
			|| inSeqFileFormat.equalsIgnoreCase("fq")){
		
		
    	 isOk=true;
	
	   }
	 
	   return isOk;
   }

   public static List<SeqInfo> checkSeq(String inSeqFile){
	    
	    List<SeqInfo> seqObjList = null;

		if(isFASTASeq(inSeqFile)){		  
			seqObjList=checkFASTASeq(inSeqFile);			
		}else if(isFASTQSeq(inSeqFile)){		  
			seqObjList=checkFASTQSeq(inSeqFile);			
		}
		
		return seqObjList;
   } 
	 
   static List<SeqInfo> checkFASTASeq(String seqFile){
		  
	    List<SeqInfo> seqObjList=new ArrayList<SeqInfo>();
	    SeqInfo perSeq;
		int seqID=0;

		try{  
			BufferedReader br;              
	        br = new BufferedReader(new FileReader(seqFile));
		    String line;
			String seqIdentifier;
			String seqLine;
			String seqName;
			String [] itemSplited;		
			
			line = br.readLine();
			if (line == null){ 
			   br.close();
			   return seqObjList;
			}
			while(line.length()==0 || line.matches("\\s*")){
			   line = br.readLine();
			   if (line == null){ 
				   br.close();
				   return seqObjList;
			   }
			}
			
			seqID=0;
			while(true){           
		       if (line == null) break;
		       line=line.trim();
			   if(line.indexOf(">")==0){			    
			     seqIdentifier=line.substring(1,line.length());			 
				 itemSplited=seqIdentifier.split("\\s+");
				 seqName=itemSplited[0].trim();
				  
				 line=br.readLine();
				 seqLine="";
				 if (line == null) break;
				 while(line.indexOf(">")<0){
			      seqLine = seqLine+line.trim();
				  line=br.readLine();
				  if (line == null) break;
				 }
				 seqLine=seqLine.replaceAll("N","n");
			  
				 perSeq=new SeqInfo();
				 perSeq.seqID=seqID;
				 perSeq.seqIdentifier=seqIdentifier;
				 perSeq.seqName=seqName;
				 perSeq.seqLength=seqLine.length();
				 perSeq.seq=seqLine;
				 seqObjList.add(perSeq);
				 perSeq=null; 
				 
				 seqID=seqID+1;
			   }
			}           		   
			br.close();
		
			line=null;
			seqIdentifier=null;	
			seqLine=null;
			seqName=null;
			itemSplited=null;
		}catch(IOException e){
	        System.out.println(e);
	    }

	    return seqObjList; 
   }
	 
   static List<SeqInfo> checkFASTQSeq(String seqFile){
	  
	    List<SeqInfo> seqObjList=new ArrayList<SeqInfo>();
	    SeqInfo perSeq;
		int seqID=0;
		
		try{    
	       	BufferedReader br;              
	        br = new BufferedReader(new FileReader(seqFile));
		    String line;
			String seqIdentifier;	
			String seqLine;
			String seqName;		
			String seqQualityLine;
			String [] itemSplited;
			seqID=0;
			outerloop:	
			while(true){ 
			   //for read header/identifier
	           line = br.readLine();			
		       if (line == null) break;
		       while(line.length()==0 || line.matches("\\s*")){
		    	 line = br.readLine();
		    	 if (line == null) break outerloop;
		       }
			   line=line.trim();
	           if(line.indexOf("@")!=0) {
			      System.out.println("Error in reading fastq line:"+line);
				  break;
			   }			   
	      	   seqIdentifier=line.substring(1,line.length());				 
			   itemSplited=seqIdentifier.split("\\s+");
			   seqName=itemSplited[0].trim();
			   
			   //for read sequence
			   line=br.readLine();			 
			   if (line == null) break;
		       while(line.length()==0 || line.matches("\\s*")){
		    	 line = br.readLine();
		    	 if (line == null) break outerloop;
		       }
			   seqLine=line.trim();			 
			   seqLine=seqLine.replaceAll("N","n");
			   
			   //for read '+' tag
			   line=br.readLine();			  
			   if (line == null) break;
		       while(line.length()==0 || line.matches("\\s*")){
		    	 line = br.readLine();
		    	 if (line == null) break outerloop;
		       }
			   line=line.trim();
			   if(line.indexOf("+")!=0) {
				  System.out.println("Error in reading fastq line:"+line);
				  break;
			   }		      
	           
			   //for read base quality encode
	           line=br.readLine();			 
			   if (line == null) break;
			   while(line.length()==0 || line.matches("\\s*")){
				  line = br.readLine();
				  if (line == null) break outerloop;
			   }
			   seqQualityLine=line.trim();

			   perSeq=new SeqInfo();
			   perSeq.seqID=seqID;
			   perSeq.seqIdentifier=seqIdentifier;
			   perSeq.seqName=seqName;
			   perSeq.seqLength=seqLine.length();
			   perSeq.seq=seqLine;
			   perSeq.seqQualityEncode=seqQualityLine;
			   seqObjList.add(perSeq);
			   perSeq=null;	
			   
			   seqID=seqID+1;	
			}           		   
			br.close();
		
	    	line=null;
	    	seqIdentifier=null;	
			seqLine=null;
			seqName=null;
			itemSplited=null;
		}catch(IOException e){
	        System.out.println(e);
	    }

	    return seqObjList; 
   }

   public static boolean checkPairEndSeq(List<SeqInfo> seqObjList, List<SeqInfo> seqObjList2){
	   
	   boolean isOK=true;
	   if(seqObjList.size()!=seqObjList2.size()){
		  System.out.println("Error: different seq num of pair-end seq");
		  return false;
	   }

	   for(int i=0;i<seqObjList.size();i++){
		 if(!seqObjList.get(i).seqName.trim().equals(seqObjList2.get(i).seqName.trim())){ 
			isOK=false;
			System.out.println("Error: inconsistent seq name between pair-end seq");
			System.out.println("Forward: "+seqObjList.get(i).seqName
					+"\t Reverse:"+seqObjList2.get(i).seqName);
			break;
		 }
	   }	
	 
	   return isOK;
   }


   public static List<SeqInfo> getSeqObj(String seqFile){
     
	   List<SeqInfo> seqObj=new ArrayList<SeqInfo>();
	   if(isFASTASeq(seqFile)){		
		 seqObj=getFASTASeqObj(seqFile);		 
	   }else if(isFASTQSeq(seqFile)){		 
		 seqObj=getFASTQSeqObj(seqFile);		 
	   }
	 
	   return seqObj;
		
   }
 
   public static List<SeqInfo> getFASTASeqObj(String seqFile){
	 
	    List<SeqInfo> seqObjList=new ArrayList<SeqInfo>();
		SeqInfo seqObj;
		int seqID=0;
		try{    
	        BufferedReader br;              
	        br = new BufferedReader(new FileReader(seqFile));
		    String line;
			String seqIdentifier;
			String seqName;
			String seqLine;
			String [] itemSplited;			
			seqID=0;
			line = br.readLine();
			if (line == null){ 
			   br.close();
			   return seqObjList;
			}
			while(line.length()==0 || line.matches("\\s*")){
			   line = br.readLine();
			   if (line == null){ 
				   br.close();
				   return seqObjList;
			   }
			}
			
			while(true){           
		       if (line == null) break;
		       line=line.trim();
			   if(line.indexOf(">")==0){			     
				 seqIdentifier=line.substring(1,line.length());			 
				 itemSplited=seqIdentifier.split("\\s+");
				 seqName=itemSplited[0].trim();	
				 
				 line=br.readLine();
				 seqLine="";
				 if (line == null) break;
				 while(line.indexOf(">")<0){
			      seqLine = seqLine+line.trim();
				  line=br.readLine();
				  if (line == null) break;
				 }
				 seqLine=seqLine.replaceAll("N","n");
	             seqObj=new SeqInfo();
	             seqObj.seqID=seqID;
	             seqObj.seqIdentifier=seqIdentifier;
	             seqObj.seqName=seqName;
				 seqObj.seqLength=seqLine.length();
	             seqObj.seq=seqLine.toUpperCase();
	             seqObjList.add(seqObj);
	             seqObj=null;
	             
	             seqID=seqID+1;
	             
			   }
				 
			}           		   
			br.close();
		
		}catch(IOException e){
	        System.out.println(e);
	    }
		
		return seqObjList;
	   
   }
	 
   public static List<SeqInfo> getFASTASeqObj(String seqFile, int start, int end){
		 
		    List<SeqInfo> seqObjList=new ArrayList<SeqInfo>();
		    SeqInfo perSeq=new SeqInfo();
			int seqID=0;
			
			try{    
		       	BufferedReader br;              
		        br = new BufferedReader(new FileReader(seqFile));
			    String line;
				String seqIdentifier="";
				String seqLine = "";
				String seqName="";			
				String [] itemSplited;
				seqID=0;
				line = br.readLine();
				if (line == null){ 
				   br.close();
				   return seqObjList;
				}
				while(line.length()==0 || line.matches("\\s*")){
				   line = br.readLine();
				   if (line == null){ 
					   br.close();
				       return seqObjList;
				   }
				}
				
				while(true){		         		
			       if (line == null) break;
			       if(seqID>=end) break;
			       line=line.trim();
				   if(line.indexOf(">")==0){			    
					  seqIdentifier=line.substring(1,line.length());				 
					  line=br.readLine();				
					  if (line == null) break;
					  seqLine="";
					  while(line.indexOf(">")<0){
				        seqLine = seqLine+line.trim();
					    line=br.readLine();
					    if (line == null) break;
					  }
					 
					  if(seqID>=start && seqID<end){					   
						 itemSplited=seqIdentifier.split("\\s+");
						 seqName=itemSplited[0].trim();					
						 seqLine=seqLine.replaceAll("N","n");
						 perSeq=new SeqInfo();				  
						 perSeq.seqIdentifier=seqIdentifier;
						 perSeq.seqID=seqID;
						 perSeq.seqName=seqName;
						 perSeq.seqLength=seqLine.length();
						 perSeq.seq=seqLine.toUpperCase();
				         seqObjList.add(perSeq);
				         perSeq=null;	
					  }
					 
					  seqID=seqID+1;
		             
				   }
				}           		   
				br.close();
			
		    	line=null;
		    	seqIdentifier=null;		
				seqLine=null;
				seqName=null;
				itemSplited=null;
			}catch(IOException e){
		        System.out.println(e);
		    }

		    return seqObjList; 
   }
	 
   public static List<SeqInfo> getFASTQSeqObj(String seqFile){
		 
		    List<SeqInfo> seqObjList=new ArrayList<SeqInfo>();
		    SeqInfo perSeq=new SeqInfo();
			int seqID=0;
			
			try{    
		       	BufferedReader br;              
		        br = new BufferedReader(new FileReader(seqFile));
			    String line;
				String seqIdentifier;	
				String seqLine;
				String seqName;		
				String seqQualityLine;
				String [] itemSplited;
				outerloop:
				while(true){  			   
				   //for read header/identifier
		           line = br.readLine();			
			       if(line == null) break;
			       while(line.length()==0 || line.matches("\\s*")){
			    	 line = br.readLine();
			    	 if (line == null) break outerloop;
			       }
				   line=line.trim();
		           if(line.indexOf("@")!=0) {
				      System.out.println("Error in reading fastq line: "+line);
					  break;
				   }			   
		      	   seqIdentifier=line.substring(1,line.length());				 
				   itemSplited=seqIdentifier.split("\\s+");
				   seqName=itemSplited[0].trim();
				   
				   //for read sequence
				   line=br.readLine();			 
				   if (line == null) break;
			       while(line.length()==0 || line.matches("\\s*")){
			    	 line = br.readLine();
			    	 if (line == null) break outerloop;
			       }
				   seqLine=line.trim();			 
				   seqLine=seqLine.replaceAll("N","n");
				   
				   //for read '+' tag
				   line=br.readLine();			  
				   if (line == null) break;
			       while(line.length()==0 || line.matches("\\s*")){
			    	 line = br.readLine();
			    	 if (line == null) break outerloop;
			       }
				   line=line.trim();
				   if(line.indexOf("+")!=0) {
					  System.out.println("Error in reading fastq line: "+line);
					  break;
				   }		      
		           
				   //for read base quality encode
		           line=br.readLine();			 
				   if (line == null) break;
				   while(line.length()==0 || line.matches("\\s*")){
					  line = br.readLine();
					  if (line == null) break outerloop;
				   }
				   seqQualityLine=line.trim();
				   
				   perSeq=new SeqInfo();
				   perSeq.seqID=seqID;
				   perSeq.seqIdentifier=seqIdentifier;
				   perSeq.seqName=seqName;
				   perSeq.seqLength=seqLine.length();
				   perSeq.seq=seqLine;
				   perSeq.seqQualityEncode=seqQualityLine;
				   seqObjList.add(perSeq);
				   perSeq=null;
				   
				   seqID=seqID+1;
			
				}           		   
				br.close();
			
		    	line=null;
		    	seqIdentifier=null;		
				seqLine=null;
				seqName=null;
				itemSplited=null;
			}catch(IOException e){
		        System.out.println(e);
		    }

		    return seqObjList; 
   }
	 
   public static List<SeqInfo> getFASTQSeqObj(String seqFile, int start, int end){
		 
		    List<SeqInfo> seqObjList=new ArrayList<SeqInfo>();
		    SeqInfo perSeq=new SeqInfo();
			int seqID=0;
			
			try{    
		       	BufferedReader br;              
		        br = new BufferedReader(new FileReader(seqFile));
			    String line;
				String seqIdentifier;		
				String seqLine;
				String seqName;		
				String seqQualityLine;
				String [] itemSplited;
				outerloop:	
				while(true){  
				   
			       if(seqID>=end) break;
			       
			       //for read header/identifier
		           line = br.readLine();			
			       if (line == null) break;
			       while(line.length()==0 || line.matches("\\s*")){
			    	 line = br.readLine();
			    	 if (line == null) break outerloop;
			       }
				   line=line.trim();
		           if(line.indexOf("@")!=0) {
				      System.out.println("Error in reading fastq line:"+line);
					  break;
				   }			   
		      	   seqIdentifier=line.substring(1,line.length());				 
				   itemSplited=seqIdentifier.split("\\s+");
				   seqName=itemSplited[0].trim();
				   
				   //for read sequence
				   line=br.readLine();			 
				   if (line == null) break;
			       while(line.length()==0 || line.matches("\\s*")){
			    	 line = br.readLine();
			    	 if (line == null) break outerloop;
			       }
				   seqLine=line.trim();			 
				   //seqLine=seqLine.replaceAll("N","n");
				   
				   //for read '+' tag
				   line=br.readLine();			  
				   if (line == null) break;
			       while(line.length()==0 || line.matches("\\s*")){
			    	 line = br.readLine();
			    	 if (line == null) break outerloop;
			       }
				   line=line.trim();
				   if(line.indexOf("+")!=0) {
					  System.out.println("Error in reading fastq line:"+line);
					  break;
				   }		      
		           
				   //for read base quality encode
		           line=br.readLine();			 
				   if (line == null) break;
				   while(line.length()==0 || line.matches("\\s*")){
					  line = br.readLine();
					  if (line == null) break outerloop;
				   }
				   seqQualityLine=line.trim();
				   
				   if(seqID>=start && seqID<end){
					   				 
					  itemSplited=seqIdentifier.split("\\s+");
					  seqName=itemSplited[0].trim();				
					  seqLine=seqLine.replaceAll("N","n");
					   
					  perSeq=new SeqInfo();
					  perSeq.seqID=seqID;
					  perSeq.seqIdentifier=seqIdentifier;
					  perSeq.seqName=seqName;
					  perSeq.seqLength=seqLine.length();
					  perSeq.seq=seqLine.trim();
					  perSeq.seqQualityEncode=seqQualityLine.trim();
					  seqObjList.add(perSeq);
					  perSeq=null;
				   }
				   
				   seqID=seqID+1;
			
				}           		   
				br.close();
			
		    	line=null;
		    	seqIdentifier=null;	
				seqLine=null;
				seqName=null;
				itemSplited=null;
			}catch(IOException e){
		        System.out.println(e);
		    }

		    return seqObjList; 
   }
 

   public static List<String> splitSeqFile(String inSeqFile, int step, String outDir){
	     
	    List<String> splitedSeqFiles=new ArrayList<String>();	
	    String timeStamp=new SimpleDateFormat("yyyyMMdd_HHmmss").format(Calendar.getInstance().getTime());
	    System.out.println("###### Spliting reads for ["+inSeqFile+"](time: "+timeStamp+").........");
	
	    if(isFASTASeq(inSeqFile)){			
		   splitedSeqFiles=splitFASTAFile(inSeqFile, step, outDir);		
		}else if(isFASTQSeq(inSeqFile)){			 
		   splitedSeqFiles=splitFASTQFile(inSeqFile, step, outDir);
		}
	     
	    timeStamp=new SimpleDateFormat("yyyyMMdd_HHmmss").format(Calendar.getInstance().getTime());
	    System.out.println("###### Reads split done for ["+inSeqFile+"](time: "+timeStamp+")");
	     
	    return splitedSeqFiles;
   }
 
   public static List<String> splitFASTAFile(String inSeqFile, int step, String outDir){
	    	   
		List<String> splitedFiles=new ArrayList<String>();
	    String outTag=inSeqFile.substring(inSeqFile.lastIndexOf("/")+1, inSeqFile.lastIndexOf("."));
		if(step<0) step=500000;
		if(outDir==null) 
		  outDir=inSeqFile.substring(0,inSeqFile.lastIndexOf("/"))+"/"
		       +inSeqFile.substring(inSeqFile.lastIndexOf("/")+1,inSeqFile.lastIndexOf("."))
		       +"_Splited";
		FileOperation.newFolder(outDir);
		int total=SeqOperation.getFASTASeqNum(inSeqFile);
		System.out.println("Total reads to be splited: "+total);    	
		try{
		  BufferedWriter writer=null;
		  int stopPos=0;		 
		  for(int f=0;f<total;f=f+step){  
	         String seqIdentifier;
	         String seq;	      
	         stopPos=Math.min(total, (f+step));
	         List<SeqInfo> seqObjList=getFASTASeqObj(inSeqFile,f,stopPos);
	         String outSeqFile=outDir+"/"+outTag+".split."+step+"_"+(f+1)+"-"+stopPos+".fna";	 
	         writer=new BufferedWriter(new FileWriter(outSeqFile));		         	
	         for(int i=0;i<seqObjList.size();i++){	
	        	seqIdentifier=seqObjList.get(i).seqIdentifier;
	 			seq=seqObjList.get(i).seq;
	 			writer.write(">"+seqIdentifier);
	 			writer.newLine();
	 			//writer.flush(); 
	 			writer.write(seq);
	 			writer.newLine();
	 			writer.flush();				
		     }
		     writer.close();
		     splitedFiles.add(outSeqFile);
		     outSeqFile=null;
		     seqObjList=null;
	         seqIdentifier=null;
	         seq=null;	    
		     System.out.println((f+1)+"-"+stopPos+" split OK!");	
		  }
		  writer=null;
		}catch(IOException e){
	        System.out.println(e);
	    }	
		  
		return  splitedFiles;
   }
 
   public static List<String> splitFASTQFile(String inSeqFile, int step, String outDir){
	    	   
		List<String> splitedFiles=new ArrayList<String>();
		String outTag=inSeqFile.substring(inSeqFile.lastIndexOf("/")+1, inSeqFile.lastIndexOf("."));
		if(step<0) step=500000;
		if(outDir==null) 
			outDir=inSeqFile.substring(0,inSeqFile.lastIndexOf("/"))+"/"
		           +inSeqFile.substring(inSeqFile.lastIndexOf("/")+1,inSeqFile.lastIndexOf("."))
		           +"_Splited";
		FileOperation.newFolder(outDir);
		int total=SeqOperation.getFASTQSeqNum(inSeqFile);
		System.out.println("Total reads to be splited: "+total);    	
		try{
		  BufferedWriter writer=null;
		  int stopPos=0;		 
		  for(int f=0;f<total;f=f+step){
	         String seqIdentifier;
	         String seq;
	         String seqQuality;
	         stopPos=Math.min(total, (f+step));
	         List<SeqInfo> seqObjList=getFASTQSeqObj(inSeqFile,f,stopPos);
	         String outSeqFile=outDir+"/"+outTag+".split."+step+"_"+(f+1)+"-"+stopPos+".fastq";	  
		     writer=new BufferedWriter(new FileWriter(outSeqFile));		         	
	         for(int i=0;i<seqObjList.size();i++){	
	        	seqIdentifier=seqObjList.get(i).seqIdentifier;
				seq=seqObjList.get(i).seq;
				seqQuality=seqObjList.get(i).seqQualityEncode;
				
				writer.write("@"+seqIdentifier);
				writer.newLine();		
				writer.write(seq);
				writer.newLine();
				writer.write("+");
				writer.newLine();
				writer.write(seqQuality);
				writer.newLine();
				writer.flush();			
		     }
		     writer.close();
		     splitedFiles.add(outSeqFile);
		     outSeqFile=null;
		     seqObjList=null;
	         seqIdentifier=null;
	         seq=null;
	         seqQuality=null;
		     System.out.println((f+1)+"-"+stopPos+" split OK!");	
		  }
		  writer=null;
		}catch(IOException e){
	        System.out.println(e);
	    }
		
		  
		return  splitedFiles;
   }
   
   
   public static String filterSeqFile(String inSeqFile, int step, int minSeqLen, String outDir){	
	     
	    String outSeqFile = null;
	    if(isFASTASeq(inSeqFile)){			
	       outSeqFile=filterFASTAFile(inSeqFile, step, minSeqLen, outDir);		
		}else if(isFASTQSeq(inSeqFile)){			 
		   outSeqFile=filterFASTQFile(inSeqFile, step, minSeqLen, outDir);
		}
	     
	    return outSeqFile;
	
   }

   public static String filterFASTAFile(String inSeqFile, int step, int minSeqLen, String outDir){		
	   	
		if(step<0) step=1000000;
		if(outDir==null) outDir=inSeqFile.substring(0,inSeqFile.lastIndexOf("/"))+"/filtered";
		FileOperation.newFolder(outDir);
		int total=SeqOperation.getFASTASeqNum(inSeqFile);
		System.out.println("Total Seq to be filtered: "+total); 
		String timeStamp=new SimpleDateFormat("yyyyMMdd_HHmmss").format(Calendar.getInstance().getTime());
		String outSeqFile = null;
		try{
		  BufferedWriter writer=null;
	      outSeqFile=outDir+"/"
		       +inSeqFile.substring(inSeqFile.lastIndexOf("/")+1,inSeqFile.lastIndexOf("."))
		       +".minSeqLen"+minSeqLen+"bp."+timeStamp+".fna";	 
	      writer=new BufferedWriter(new FileWriter(outSeqFile));
		  int stopPos=0;		 
		  for(int f=0;f<total;f=f+step){  
	         String seqIdentifier;
	         String seq;	      
	         stopPos=Math.min(total, (f+step));	
	         List<SeqInfo> seqObjList=getFASTASeqObj(inSeqFile,f,stopPos);	
	         for(int i=0;i<seqObjList.size();i++){	
	        	seqIdentifier=seqObjList.get(i).seqIdentifier;
	 			seq=seqObjList.get(i).seq;
	 			if(seq.length()>=minSeqLen) {
	 			  writer.write(">"+seqIdentifier);
	 		 	  writer.newLine();
	 			  //writer.flush(); 
	 			  writer.write(seq);
	 			  writer.newLine();
	 			  writer.flush();	
	            }
		     }
	         seqObjList=null;
	         seqIdentifier=null;
	         seq=null;	     
		     System.out.println((f+1)+"-"+stopPos+" filter OK!");	
		  }
		  writer.close();
		  writer=null;
		  outSeqFile=null;
		}catch(IOException e){
	      System.out.println(e);
	    }
		
		timeStamp=new SimpleDateFormat("yyyyMMdd_HHmmss").format(Calendar.getInstance().getTime());
		System.out.println("Seq filter done! [time:"+timeStamp+"]");	
		timeStamp=null;		
		
		return outSeqFile;

   }

   public static String filterFASTQFile(String inSeqFile, int step, int minSeqLen, String outDir){	    	   
	    
		if(step<0) step=1000000;
		if(outDir==null) outDir=inSeqFile.substring(0,inSeqFile.lastIndexOf("/"))+"/filtered";
		FileOperation.newFolder(outDir);
		int total=SeqOperation.getFASTQSeqNum(inSeqFile);
		System.out.println("Total Seq to be filtered: "+total);    	
		String timeStamp=new SimpleDateFormat("yyyyMMdd_HHmmss").format(Calendar.getInstance().getTime());
		String outSeqFile=null;
		try{
		  BufferedWriter writer=null;
	      outSeqFile=outDir+"/"
			       +inSeqFile.substring(inSeqFile.lastIndexOf("/")+1,inSeqFile.lastIndexOf("."))
			       +".minSeqLen"+minSeqLen+"bp."+timeStamp+".fastq";	
	      writer=new BufferedWriter(new FileWriter(outSeqFile));
		  int stopPos=0;	
		 
		  for(int f=0;f<total;f=f+step){
	         String seqIdentifier;
	         String seq;
	         String seqQuality;
	         stopPos=Math.min(total, (f+step));
	         List<SeqInfo> seqObjList=getFASTQSeqObj(inSeqFile,f,stopPos);	
	         for(int i=0;i<seqObjList.size();i++){	
	        	seqIdentifier=seqObjList.get(i).seqIdentifier;
				seq=seqObjList.get(i).seq;
				seqQuality=seqObjList.get(i).seqQualityEncode;
				if(seq.length()>=minSeqLen) {
				  writer.write("@"+seqIdentifier);
				  writer.newLine();		
				  writer.write(seq);
				  writer.newLine();
				  writer.write("+");
				  writer.newLine();
				  writer.write(seqQuality);
				  writer.newLine();
				  writer.flush();	
				}
		     }
	         seqObjList=null;
	         seqIdentifier=null;
	         seq=null;
	         seqQuality=null;
		     System.out.println((f+1)+"-"+stopPos+" filter OK!");	
		  }
		  writer.close();
		  writer=null;
		  outSeqFile=null;
		}catch(IOException e){
	      System.out.println(e);
	    }
		
		timeStamp=new SimpleDateFormat("yyyyMMdd_HHmmss").format(Calendar.getInstance().getTime());
		System.out.println("Seq filter done! [time:"+timeStamp+"]");	
		timeStamp=null;
		
		return outSeqFile;

   }   
   
   public static String[] filterSeqFile(String inSeqFile, String inSeqFile2, int step, int minSeqLen, String outDir){	
	     
	     String[] outSeqFile12=new String[2];
	     
	     if(isFASTASeq(inSeqFile)){			
	    	 outSeqFile12=filterFASTAFile(inSeqFile, inSeqFile2, step, minSeqLen, outDir);		
		 }else if(isFASTQSeq(inSeqFile)){			 
			filterFASTQFile(inSeqFile, inSeqFile2, step, minSeqLen, outDir);
		 }
	     
	     return outSeqFile12;
   }

   public static String[] filterFASTAFile(String inSeqFile, String inSeqFile2, int step, int minSeqLen, String outDir){		
	   	
		if(step<0) step=1000000;
		if(outDir==null) outDir=inSeqFile.substring(0,inSeqFile.lastIndexOf("/"))+"/filtered";
		FileOperation.newFolder(outDir);
		int total=SeqOperation.getFASTASeqNum(inSeqFile);
		System.out.println("Total Seq to be filtered: "+total); 
		String timeStamp=new SimpleDateFormat("yyyyMMdd_HHmmss").format(Calendar.getInstance().getTime());
		String[] outSeqFile12=new String[2];
		try{
		  BufferedWriter writer=null;
	      String outSeqFile=outDir+"/"
		       +inSeqFile.substring(inSeqFile.lastIndexOf("/")+1,inSeqFile.lastIndexOf("."))
		       +".minSeqLen"+minSeqLen+"bp."+timeStamp+".fna";	 
	      writer=new BufferedWriter(new FileWriter(outSeqFile));
	      
		  BufferedWriter writer2=null;
	      String outSeqFile2=outDir+"/"
		       +inSeqFile2.substring(inSeqFile2.lastIndexOf("/")+1,inSeqFile2.lastIndexOf("."))
		       +".minSeqLen"+minSeqLen+"bp."+timeStamp+".fna";	 
	      writer2=new BufferedWriter(new FileWriter(outSeqFile2));
	      
	      outSeqFile12[0]=outSeqFile;
	      outSeqFile12[1]=outSeqFile2;
	      
		  int stopPos=0;		 
		  for(int f=0;f<total;f=f+step){  
	         String seqIdentifier;
	         String seq;
	         String seqIdentifier2;
	         String seq2;	 
	         stopPos=Math.min(total, (f+step));	
	         List<SeqInfo> seqObjList=getFASTASeqObj(inSeqFile,f,stopPos);
	         List<SeqInfo> seqObjList2=getFASTASeqObj(inSeqFile2,f,stopPos);
	         for(int i=0;i<seqObjList.size();i++){	
	        	
	        	seqIdentifier=seqObjList.get(i).seqIdentifier;
	 			seq=seqObjList.get(i).seq;
	 			
	        	seqIdentifier2=seqObjList2.get(i).seqIdentifier;
	 			seq2=seqObjList2.get(i).seq;
	 			
	 			if(seq.length()>=minSeqLen || seq2.length()>=minSeqLen) {
	 			  //forward
	 			  writer.write(">"+seqIdentifier);
	 		 	  writer.newLine();
	 			  //writer.flush(); 
	 			  writer.write(seq);
	 			  writer.newLine();
	 			  writer.flush();	
	 			  
	 			  //reverse
	 			  writer2.write(">"+seqIdentifier2);
	 		 	  writer2.newLine();
	 			  //writer2.flush(); 
	 			  writer2.write(seq2);
	 			  writer2.newLine();
	 			  writer2.flush();	

	            }
	 			
		     }
	 		 seqObjList=null;
	 		 seqObjList2=null;
	         seqIdentifier=null;
	         seq=null;
	         seqIdentifier2=null;
	         seq2=null;	
		     System.out.println((f+1)+"-"+stopPos+" filter OK!");	
		  }
		  writer.close();		  
		  outSeqFile=null;
		  writer2.close();		  
		  outSeqFile2=null;
		}catch(IOException e){
	      System.out.println(e);
	    }

		timeStamp=new SimpleDateFormat("yyyyMMdd_HHmmss").format(Calendar.getInstance().getTime());
		System.out.println("Seq filter done! [time:"+timeStamp+"]");	
		timeStamp=null;
		
		return outSeqFile12;

   }

   public static String[] filterFASTQFile(String inSeqFile, String inSeqFile2, int step, int minSeqLen, String outDir){	    	   
	    
		if(step<0) step=1000000;
		if(outDir==null) outDir=inSeqFile.substring(0,inSeqFile.lastIndexOf("/"))+"/filtered";
		FileOperation.newFolder(outDir);
		int total=SeqOperation.getFASTQSeqNum(inSeqFile);
		System.out.println("Total Seq to be filtered: "+total);    	
		String timeStamp=new SimpleDateFormat("yyyyMMdd_HHmmss").format(Calendar.getInstance().getTime());
		String[] outSeqFile12=new String[2];
		
		try{
		  BufferedWriter writer=null;
	      String outSeqFile=outDir+"/"
			       +inSeqFile.substring(inSeqFile.lastIndexOf("/")+1,inSeqFile.lastIndexOf("."))
			       +".minSeqLen"+minSeqLen+"bp."+timeStamp+".fastq";	
	      writer=new BufferedWriter(new FileWriter(outSeqFile));
	      
		  BufferedWriter writer2=null;
	      String outSeqFile2=outDir+"/"
			       +inSeqFile2.substring(inSeqFile2.lastIndexOf("/")+1,inSeqFile2.lastIndexOf("."))
			       +".minSeqLen"+minSeqLen+"bp."+timeStamp+".fastq";	
	      writer2=new BufferedWriter(new FileWriter(outSeqFile2));
	      
	      outSeqFile12[0]=outSeqFile;
	      outSeqFile12[1]=outSeqFile2;
	      
		  int stopPos=0;		 
		  for(int f=0;f<total;f=f+step){
	         String seqIdentifier;
	         String seq;
	         String seqQuality;
	         String seqIdentifier2;
	         String seq2;
	         String seqQuality2;
	         stopPos=Math.min(total, (f+step));
	         List<SeqInfo> seqObjList=getFASTQSeqObj(inSeqFile,f,stopPos);	
	         List<SeqInfo> seqObjList2=getFASTQSeqObj(inSeqFile,f,stopPos);	
	         for(int i=0;i<seqObjList.size();i++){	
	        	
	        	seqIdentifier=seqObjList.get(i).seqIdentifier;
				seq=seqObjList.get(i).seq;
				seqQuality=seqObjList.get(i).seqQualityEncode;
				
	        	seqIdentifier2=seqObjList2.get(i).seqIdentifier;
				seq2=seqObjList2.get(i).seq;
				seqQuality2=seqObjList2.get(i).seqQualityEncode;
				
				if(seq.length()>=minSeqLen || seq2.length()>=minSeqLen) {
				  //forward
				  writer.write("@"+seqIdentifier);
				  writer.newLine();		
				  writer.write(seq);
				  writer.newLine();
				  writer.write("+");
				  writer.newLine();
				  writer.write(seqQuality);
				  writer.newLine();
				  writer.flush();	
				  
				  //reverse
				  writer2.write("@"+seqIdentifier2);
				  writer2.newLine();		
				  writer2.write(seq2);
				  writer2.newLine();
				  writer2.write("+");
				  writer2.newLine();
				  writer2.write(seqQuality2);
				  writer2.newLine();
				  writer2.flush();	

				}				

		     }
		     seqIdentifier=null;
			 seq=null;
			 seqQuality=null;
					
		     seqIdentifier2=null;
			 seq2=null;
			 seqQuality2=null;

	         seqObjList=null;
	         seqObjList2=null;
		     System.out.println((f+1)+"-"+stopPos+" filter OK!");	
		  }
		  writer.close();	
		  writer=null;
		  writer2.close();	
		  writer2=null;
		  outSeqFile=null;
		  outSeqFile2=null;

		}catch(IOException e){
	      System.out.println(e);
	    }
		
		timeStamp=new SimpleDateFormat("yyyyMMdd_HHmmss").format(Calendar.getInstance().getTime());
		System.out.println("Seq filter done! [time:"+timeStamp+"]");	
		timeStamp=null;
		
		return outSeqFile12;

   }

   
   public static String convertFASTQ2FASTA(String inSeqFile){
		  
	    List<SeqInfo> seqObjList;	
	    int total=SeqOperation.getFASTQSeqNum(inSeqFile);
	    String outSeqFile=inSeqFile.substring(0,inSeqFile.lastIndexOf("."))+".fasta";	
		try{
		  BufferedWriter writer=null;		   
		  writer=new BufferedWriter(new FileWriter(outSeqFile));	
		  int stopPos=0;
		  for(int f=0;f<total;f=f+splitStep){	
	         String seqIdentifier;
	         String seq;	
	         stopPos=Math.min(total, (f+splitStep));
	         seqObjList=getFASTQSeqObj(inSeqFile,f,stopPos);	
	         for(int i=0;i<seqObjList.size();i++){	
	        	seqIdentifier=seqObjList.get(i).seqIdentifier;
				seq=seqObjList.get(i).seq;					
				writer.write(">"+seqIdentifier);
				writer.newLine();		
				writer.write(seq);
				writer.newLine();

				writer.flush();			
		     }
		  }
		  writer.close();
		}catch(IOException e){
	      System.out.println(e);
	    }
		seqObjList=null;
		
		return outSeqFile;

   }  


   public static String combineSeqFile(List<String> splitedFiles, String outFile){
	 
	   if(isFASTASeq(splitedFiles)){		
		  outFile=combineFASTAFile(splitedFiles, outFile);		 
	   }else if(isFASTQSeq(splitedFiles)){		 
		  outFile=combineFASTQFile(splitedFiles, outFile);		 
	   }	
	 
	   return outFile;
   }
 
   public static String combineFASTAFile(List<String> splitedFiles, String outFile){
	    
	    if(splitedFiles==null || splitedFiles.size()==0) return null;  
	   
		if(outFile==null){
			String fastaFile=splitedFiles.get(0);
			String outDir=fastaFile.substring(0,fastaFile.lastIndexOf("/"))+"/combined";
			outFile=outDir+"/"
			        +fastaFile.substring(fastaFile.lastIndexOf("/"),fastaFile.lastIndexOf("."))
			        +"_combined.fna";
			FileOperation.newFolder(outDir);
		}		
		//System.out.println("Combining.................. ");	
		try{
		  BufferedWriter writer=null;
		  writer=new BufferedWriter(new FileWriter(outFile));
		  String seqIdentifier;
		  String seq;	
		  int num=0;
		  for(String file:splitedFiles){	
			 if(file==null || !new File(file).exists()) continue;
			 List<SeqInfo> seqObjList=getSeqObj(file);	 
			 num=num+seqObjList.size();
		     for(int i=0;i<seqObjList.size();i++){		
		    	seqIdentifier=seqObjList.get(i).seqIdentifier;
				seq=seqObjList.get(i).seq;
				writer.write(">"+seqIdentifier);
				writer.newLine();			
				writer.write(seq);
				writer.newLine();
				writer.flush();			
		     }	     
		     seqObjList=null;
		  }
		  writer.close();
		  writer=null;
		  System.out.println("Totally combined number: "+num);

		}catch(IOException e){
	      System.out.println(e);
	    }
			  
		return  outFile;
   }
 
   public static String combineFASTQFile(List<String> splitedFiles, String outFile){
		
	    if(splitedFiles==null || splitedFiles.size()==0) return null;
		if(outFile==null) {
		   String fastqFile=splitedFiles.get(0);
		   String outDir=fastqFile.substring(0,fastqFile.lastIndexOf("/"))+"/combined";
		   outFile=outDir+"/"
				        +fastqFile.substring(fastqFile.lastIndexOf("/"),fastqFile.lastIndexOf("."))
				        +"_combined.fastq";
		   FileOperation.newFolder(outDir);
		}
		
	   //System.out.println("Combining.................. ");

		try{			
		  BufferedWriter writer=null;
		  writer=new BufferedWriter(new FileWriter(outFile));			  
		  String seqIdentifier;
		  String seq;	
		  String seqQuality;  			
		  int num=0;		  
		  for(String file:splitedFiles ){	
			 if(file==null || !new File(file).exists()) continue;
			 List<SeqInfo> seqObjList=getSeqObj(file);	 
			 num=num+seqObjList.size();
		     for(int i=0;i<seqObjList.size();i++){		
		    	seqIdentifier=seqObjList.get(i).seqIdentifier;
				seq=seqObjList.get(i).seq;
				seqQuality=seqObjList.get(i).seqQualityEncode;
				
				writer.write("@"+seqIdentifier);
				writer.newLine();		
				writer.write(seq);
				writer.newLine();
				writer.write("+");
				writer.newLine();
				writer.write(seqQuality);
				writer.newLine();
				writer.flush();		
				writer.flush();			
		     }	     
		     seqObjList=null;
		  }
		  writer.close();
		  writer=null;
		  System.out.println("Totally combined number: "+num);

		}catch(IOException e){
	      System.out.println(e);
	    }
			  
		return  outFile;
   }	 
 
   public static int getSeqNum(List<SeqInfo> seqObjList){
  
      return seqObjList.size();	
	
   }
 
   public static int getSeqNum(String seqFile){
	 
	   int seqNum=0;

	   if(isFASTASeq(seqFile)){		
		 seqNum=getFASTASeqNum(seqFile);		 
	   }else if(isFASTQSeq(seqFile)){		 
		 seqNum=getFASTQSeqNum(seqFile);		 
	   }
	 
	   return seqNum;
		
   }
 
   public static int getFASTASeqNum(String seqFile){
  
      int seqNum=0;
	  try{    
		BufferedReader br;              
        br = new BufferedReader(new FileReader(seqFile));
	    String line;		
		seqNum=0;						
		while(true){
           line = br.readLine();
	       if (line == null) break;
		   if(line.trim().indexOf(">")==0){
		     seqNum=seqNum+1;			
		   }			 
		}           		   
		br.close();		
	  }catch(IOException e){
        System.out.println(e);
      }   
	
      return seqNum;	  
   }
 
   public static int getFASTQSeqNum(String seqFile){
	  
	    int seqNum=0;
		try{    
			BufferedReader br;              
	        br = new BufferedReader(new FileReader(seqFile));
		    String line;		
			seqNum=0;						
			while(true){
	           line = br.readLine(); // seq name line
		       if (line == null) break;
			   if(line.trim().indexOf("@")==0){
			     seqNum=seqNum+1;
			     line = br.readLine(); // seq line
				 line = br.readLine(); // tag '+' line
				 line = br.readLine(); // seq quality line		
			   }			 
			}           		   
			br.close();		
		}catch(IOException e){
	        System.out.println(e);
	    }   
		
	    return seqNum;	  
   }
 
   public static boolean extratSeqName(String seqFile, String outFile){
	    
	    boolean isOK=true;
	
	   	try{
		  if(outFile==null) outFile=seqFile+".seqName";
		  if(isFASTASeq(seqFile)){			  
			 isOK=getFASTASeqName(seqFile,outFile);			 
		  }else if(isFASTQSeq(seqFile)){			  
			 isOK=getFASTQSeqName(seqFile,outFile);			 
		  }		  
			
		}catch(Exception e){
		  isOK=false;
		  System.out.println(e);
		}
		
	   	return isOK;
	
   }
 
   public static boolean getFASTASeqName(String seqFile, String outFile){
        boolean isOK=true;
		try{    
			BufferedReader br;              
	        br = new BufferedReader(new FileReader(seqFile));
		    String line;		
		    BufferedWriter writer=null;
		    writer=new BufferedWriter(new FileWriter(outFile));
			String seqName;
			String [] itemSplited;
			while(true){
	           line = br.readLine();
		       if (line == null) break;
			   if(line.trim().indexOf(">")==0){
				 itemSplited=line.split("\\s+");
				 seqName=itemSplited[0].trim();
				 seqName=seqName.substring(line.indexOf(">")+1,seqName.length());
				 writer.write(seqName);
				 writer.newLine();
				 writer.flush();
			   }			 
			}           		   
			br.close();
			writer.close();
			writer=null;
		}catch(IOException e){
			isOK=false;
	        System.out.println(e);
	    }   
		
	    return isOK;
   }
	 
   public static boolean getFASTQSeqName(String seqFile,String outFile){
	        boolean isOK=true;
			try{    
				BufferedReader br;              
		        br = new BufferedReader(new FileReader(seqFile));
			    String line;		
			    BufferedWriter writer=null;
			    writer=new BufferedWriter(new FileWriter(outFile));
				String seqName;
				String [] itemSplited;
			
			    while(true){		         
			       line = br.readLine(); // seq name line
			       if (line == null) break;
			       line=line.trim();
				   if(line.indexOf("@")==0){
					 itemSplited=line.split("\\s+");
					 seqName=itemSplited[0].trim();
					 seqName=seqName.substring(seqName.indexOf("@")+1,seqName.length());
					 writer.write(seqName);
					 writer.newLine();
					 writer.flush();
					 line = br.readLine(); // seq line
					 line = br.readLine(); // tag '+' line
					 line = br.readLine(); // seq quality line					
				   }			 
				}           		   
				br.close();
				writer.close();
				writer=null;
			}catch(IOException e){
				isOK=false;
		        System.out.println(e);
		    } 
			
			return isOK;
		
    }
 
  public static List<String> getRowName(String nameFile){
  
    List<String> seqName=new ArrayList<String>();
	List<ArrayList <String>> nameList=FileOperation.getMatrixFromFile(nameFile);	
	String name="";
	for(int i=nameStartRowIdx; i<nameList.size();i++){
	  name=nameList.get(i).get(nameColIdx);
	  if(!seqName.contains(name)){
	    seqName.add(name);	   
	  }else{	    
	    //System.out.println("Repeated Seq: "+nameList.get(i));
	  }
	}
	nameList=null;
	System.out.println("Got Unique names: "+seqName.size() +" From ["+nameFile+"]");
	
    return seqName;
 }
 
 public static List<String> getRowName(String nameFile, int nameColIdx){
	  
	    List<String> seqName=new ArrayList<String>();
		List<ArrayList <String>> nameList=FileOperation.getMatrixFromFile(nameFile);	
		String name="";
		for(int i=nameStartRowIdx; i<nameList.size();i++){
		  name=nameList.get(i).get(nameColIdx);
		  if(!seqName.contains(name)){
		    seqName.add(name);	   
		  }else{	    
		    //System.out.println("Repeated name: "+nameList.get(i));
		  }
		}
		
		nameList=null;
		System.out.println("Got Unique names: "+seqName.size() +" From ["+nameFile+"]");
		
	    return seqName;
 }
 
 
 static List<String> getRowName(List <String> nameList,int start, int end){
	  
	    List<String> names=new ArrayList<String>();			
		String name="";
		end=Math.min(end, nameList.size());
		for(int i=start; i<end;i++){
		  name=nameList.get(i);
		  if(!names.contains(name)){
			 names.add(name);	   
		  }else{	    
		     System.out.println("Repeated name: "+nameList.get(i));
		  }
		}
		
	    return names;
 }
 
 
 public static List<SeqInfo> excludeSeq(String seqFile,List<String> excludedSeqNameList){
	    
	    List<SeqInfo> seqObjList = null;

	   	try{
			
	   	  if(isFASTASeq(seqFile)){
	   		 seqObjList=getFASTASeqObj(seqFile);
		  }else if(isFASTQSeq(seqFile)){
			 seqObjList=getFASTQSeqObj(seqFile);
		  }
		  List<String> seqNameList=new ArrayList<String>();
		  for(SeqInfo seqObj:seqObjList){
			seqNameList.add(seqObj.seqName);
		  }		
	
		  int i=0;
		  for(String name: excludedSeqNameList){
			 i=seqNameList.indexOf(name);
			 if(i>=0){
			  seqObjList.remove(i);
			  seqNameList.remove(i);
			 }
		  }
	
		  seqNameList=null;
			
		}catch(Exception e){
		  System.out.println(e);
		}
		
		return seqObjList;
 }
 
 public static List<SeqInfo> excludeSeq(String seqFile,String excludedSeqNameFile,
		 int excludedSeqNameColIdx){
	    
	    List<SeqInfo> seqObjList = null;

	   	try{
			
	   	  if(isFASTASeq(seqFile)){
			 seqObjList=getFASTASeqObj(seqFile);
		  }else if(isFASTQSeq(seqFile)){
			 seqObjList=getFASTQSeqObj(seqFile);
		  }
		  List<String> seqNameList=new ArrayList<String>();
		  for(SeqInfo seqObj:seqObjList){
			 seqNameList.add(seqObj.seqName);
		  }		
		  List<String> excludedSeqNameList=getRowName(excludedSeqNameFile,excludedSeqNameColIdx);
		  int i=0;
		  for(String name: excludedSeqNameList){
			 i=seqNameList.indexOf(name);
			 if(i>=0){
				seqObjList.remove(i);
				seqNameList.remove(i);
			 }
		  }
		
		  seqNameList=null;
		  excludedSeqNameList=null;
			
		}catch(Exception e){
		  System.out.println(e);
		}
		
		return seqObjList;
  }
 
  public static List<SeqInfo> excludeSeq(List<SeqInfo> seqObjList, List<String> excludedSeqNameList){

	   	try{
	
		  List<String> seqNameList=new ArrayList<String>();
		  for(SeqInfo seqObj:seqObjList){
			seqNameList.add(seqObj.seqName);
		  }		
	
		  int i=0;
		  for(String name: excludedSeqNameList){
			 i=seqNameList.indexOf(name);
			 if(i>=0){
			  seqObjList.remove(i);
			  seqNameList.remove(i);
			 }
		  }
	
		  seqNameList=null;
			
		}catch(Exception e){
		  System.out.println(e);
		}
		
		return seqObjList;
  }
 
  public static List<SeqInfo> extractSubSeq(String seqFile,String subSeqNameFile,int seqNameColIdx){
	    
	    List<String> totalSubSeqNames=getRowName(subSeqNameFile,seqNameColIdx);
	    int subSeqNum=totalSubSeqNames.size();
	 
	    List<SeqInfo> subSeqObjList = new ArrayList<SeqInfo>();
	    SeqInfo seqInfo;
	    int total;
	    if(isFASTASeq(seqFile))
	      total=SeqOperation.getFASTASeqNum(seqFile);
	    else if(isFASTQSeq(seqFile))
		  total=SeqOperation.getFASTQSeqNum(seqFile);
		else
		  return subSeqObjList;
	    
	    int subSeqStep=1000;
	    String seqName;
	    String seq;
		int stopPos=0;	    
	    int start;
	    int end;
	    List<String> subSeqNameList;
	   	try{
	  	    for(int f=0;f<total;f=f+splitStep){
		  		stopPos=Math.min(total, (f+splitStep));
		   		List<SeqInfo> seqObjList = null;
				if(isFASTASeq(seqFile)){
				   seqObjList=getFASTASeqObj(seqFile,f,stopPos);
				}else if(isFASTQSeq(seqFile)){
				   seqObjList=getFASTQSeqObj(seqFile,f,stopPos);
				}
				List<String> seqNameList=new ArrayList<String>();
				for(SeqInfo seqObj:seqObjList){
				  seqNameList.add(seqObj.seqName);
				}	
				
				for(start=0;start<subSeqNum; start=start+subSeqStep){			
				  end=Math.min(subSeqNum, (start+subSeqStep));
				  subSeqNameList=getRowName(totalSubSeqNames,start,end);
				  int i=0;
				  for(String name: subSeqNameList){
				    i=seqNameList.indexOf(name);
					if(i>=0){
						seqInfo=new SeqInfo();
						seqName=seqObjList.get(i).seqName;
						seq=seqObjList.get(i).seq;
						seqInfo.seqIdentifier=seqObjList.get(i).seqIdentifier;
						seqInfo.seqName=seqName;
						seqInfo.seq=seq;
						seqInfo.seqQualityEncode=seqObjList.get(i).seqQualityEncode;
						subSeqObjList.add(seqInfo);
						seqInfo=null;
					}
				  }
				  System.out.println((f+1)+"."+(start+1)+"-"+end+" OK!");
				  subSeqNameList=null;
				}
				seqObjList=null;
				seqNameList=null;
		    }
		}catch(Exception e){
		  System.out.println(e);
		}
		
		return subSeqObjList;
  } 
 
  public static void extractSubSeq(String seqFile,String subSeqNameFile,int seqNameColIdx,String outSeqFile){
	    
	    List<String> totalSubSeqNames=getRowName(subSeqNameFile,seqNameColIdx);
	    int subSeqNum=totalSubSeqNames.size();

	    int total;
	    if(isFASTASeq(seqFile))
	      total=SeqOperation.getFASTASeqNum(seqFile);
	    else if(isFASTQSeq(seqFile))
		  total=SeqOperation.getFASTQSeqNum(seqFile);
		else
		  return ;
		
		int stopPos=0;
	    int step=10000;
	    int start;
	    int end;
	    List<String> subSeqNameList;
	   	try{
		    String seqIdentifier;
		    String seq;
			String seqQuality="";
			BufferedWriter writer=null;
			writer=new BufferedWriter(new FileWriter(outSeqFile));
			
	   		for(int f=0;f<total;f=f+splitStep){
		  		stopPos=Math.min(total, (f+splitStep));
		   		List<SeqInfo> seqObjList = null;
				if(isFASTASeq(seqFile)){
				   seqObjList=getFASTASeqObj(seqFile,f,stopPos);
				}else if(isFASTQSeq(seqFile)){
				   seqObjList=getFASTQSeqObj(seqFile,f,stopPos);
				}
				List<String> seqNameList=new ArrayList<String>();
				for(SeqInfo seqObj:seqObjList){
				  seqNameList.add(seqObj.seqName);
				}	
				
				for(start=0;start<subSeqNum; start=start+step){				
				  end=Math.min(subSeqNum, (start+step));
				  subSeqNameList=getRowName(totalSubSeqNames,start,end);
			
				  if(isFASTQSeq(outSeqFile)){
						int i=0;
						for(String name: subSeqNameList){
						    i=seqNameList.indexOf(name);
							if(i>=0){
								seqIdentifier=seqObjList.get(i).seqIdentifier;
								seq=seqObjList.get(i).seq;
								seqQuality=seqObjList.get(i).seqQualityEncode;
								
								writer.write("@"+seqIdentifier);
								writer.newLine();
								writer.write(seq);
								writer.newLine();
								writer.write("+");
								writer.newLine();
								writer.write(seqQuality);
								writer.newLine();
								writer.flush();
							}
						}
				  }else{ // fasta format
						int i=0;
						for(String name: subSeqNameList){
						    i=seqNameList.indexOf(name);
							if(i>=0){
								seqIdentifier=seqObjList.get(i).seqIdentifier;
								seq=seqObjList.get(i).seq;
								writer.write(">"+seqIdentifier);
								writer.newLine();
								writer.write(seq);
								writer.newLine();
								writer.flush();
							}
						}
				  } 
				  
				  System.out.println((f+1)+"."+(start+1)+"-"+end+" OK!");
				  subSeqNameList=null;
				}
				seqObjList=null;
				seqNameList=null;
		    }
	  	    writer.close();
		}catch(Exception e){
		  System.out.println(e);
		}

  } 
 

  public static int getSubSeq(String seqFile,String subSeqNameFile,int seqNameColIdx,String outSeqFile){
	 
	 int num=0;
	 String format=FileOperation.getFileFormat(seqFile);
	 if(outSeqFile==null) outSeqFile=subSeqNameFile.substring(0,subSeqNameFile.lastIndexOf("."))+"."+format;
		
	 try{
		List<SeqInfo> seqObjList = null;
		if(isFASTASeq(seqFile)){
		   seqObjList=getFASTASeqObj(seqFile);
		}else if(isFASTQSeq(seqFile)){
		   seqObjList=getFASTQSeqObj(seqFile);
		}
		List<String> subSeqNameList=getRowName(subSeqNameFile,seqNameColIdx);
		num=getSubSeq(seqObjList,subSeqNameList,outSeqFile);
		seqObjList=null;
		subSeqNameList=null;
			
	 }catch(Exception e){
		  System.out.println(e);
	 }
	 
	 return num;
	
  } 
 
  public static int getSubSeq(List<SeqInfo> seqObjList,List<String> subSeqNameList,String outSeqFile){	
	
	 int num=0;
	 try{		
		    
		//if(seqObjList==null || seqObjList.size()==0) return num;
		//if(subSeqNameList==null || subSeqNameList.size()==0) return num;
		    
		List<String> seqNameList=new ArrayList<String>();
		for(SeqInfo seqObj:seqObjList){
		  seqNameList.add(seqObj.seqName);
		}		
			
		String seqIdentifier="";
		String seq="";
		String seqQuality="";
		BufferedWriter writer=null;
		writer=new BufferedWriter(new FileWriter(outSeqFile));
		if(isFASTQSeq(outSeqFile)){
			int i=0;
			for(String name: subSeqNameList){
			    i=seqNameList.indexOf(name);
				if(i>=0){
					seqIdentifier=seqObjList.get(i).seqIdentifier;
					seq=seqObjList.get(i).seq;
					seqQuality=seqObjList.get(i).seqQualityEncode;
						
					writer.write("@"+seqIdentifier);
					writer.newLine();
					writer.write(seq);
					writer.newLine();
					writer.write("+");
					writer.newLine();
					writer.write(seqQuality);
					writer.newLine();
					writer.flush();
					
					num++;
				}
			}
		}else{ // fasta format
			int i=0;
			for(String name: subSeqNameList){
			    i=seqNameList.indexOf(name);
				if(i>=0){
					seqIdentifier=seqObjList.get(i).seqIdentifier;
					seq=seqObjList.get(i).seq;
					writer.write(">"+seqIdentifier);
					writer.newLine();
					writer.write(seq);
					writer.newLine();
					writer.flush();
						
					num++;
				}
			}
		} 
			
		seqNameList=null;
		writer.close();
		writer=null;
			
     }catch(Exception e){
		  System.out.println(e);
	 }
	 
	 return num;
  }
 
  public static boolean isFASTAFormat(String seqFile){
	 boolean isOK=false;
	 String seqFormat=FileOperation.getFileFormat(seqFile);
	 if(seqFormat.equalsIgnoreCase("fasta") 
				|| seqFormat.equalsIgnoreCase("fna") 
				|| seqFormat.equalsIgnoreCase("fa")){
	    	isOK=true;
	 }
	 
	 return isOK;
  }
 
  public static boolean isFASTQFormat(String seqFile){
	 boolean isOK=false;
	 String seqFormat=FileOperation.getFileFormat(seqFile);
	 if(seqFormat.equalsIgnoreCase("fastq") 
				|| seqFormat.equalsIgnoreCase("fnq") 
				|| seqFormat.equalsIgnoreCase("fq")){
	    	isOK=true;
	 }
	 
	 return isOK;
  }

 
  public static void saveSeqObj(List<SeqInfo> seqObjList, String outSeqFile){	
	 if(isFASTAFormat(outSeqFile)){		
		saveSeqObjAsFASTA(seqObjList,outSeqFile);		 
	 }else if(isFASTQFormat(outSeqFile)){		 
		saveSeqObjAsFASTQ(seqObjList,outSeqFile);		 
	 }else {
		System.out.println("Warning: failed to save to ["
	          +outSeqFile+"] because seq format is not clear, please check your seq file format!!"); 
	 }		
  }
 
  public static void saveSeqObjAsFASTA(List<SeqInfo> seqObjList, String outSeqFile){
	 
	 try{
		String seqIdentifier="";
		String seq="";
		BufferedWriter writer=null;
		if(seqObjList!=null && outSeqFile!=null){
		  writer=new BufferedWriter(new FileWriter(outSeqFile));			
		  for(SeqInfo seqInfo : seqObjList){		
			seqIdentifier=seqInfo.seqIdentifier;
			seq=seqInfo.seq;
			writer.write(">"+seqIdentifier);
			writer.newLine();	
			writer.write(seq);
			writer.newLine();
			writer.flush();				
		  }
		  writer.close();
		  writer=null;
		}else{
		  System.err.println("null SeqInfo or outSeqFile!"); 
		}	
	 }catch(Exception e){
		  System.out.println(e);
	 }	 
  }
 
  public static void saveSeqObjAsFASTA(SeqInfo seqInfo, String outSeqFile){
	 
	 try{
		String seqName="";
		String seq="";
		BufferedWriter writer=null;
		if(seqInfo!=null && outSeqFile!=null){
		    writer=new BufferedWriter(new FileWriter(outSeqFile));
		
		    seqName=seqInfo.seqName;
			seq=seqInfo.seq;
			writer.write(">"+seqName+" "+seq.length());
			writer.newLine();	
			writer.write(seq);
			writer.newLine();
			writer.flush();				
		
		    writer.close();
		    writer=null;
	   }else{
		  System.err.println("Empty SeqInfo or outSeqFile!"); 
	   }	
	 }catch(Exception e){
		  System.out.println(e);
	 }	 
  }
 
  public static void saveSeqObjAsFASTQ(List<SeqInfo> seqObjList, String outSeqFile){
	 
	 try{
		String seqIdentifier="";
		String seq="";
		String seqQuality="";
		BufferedWriter writer=null;
		if(seqObjList!=null && outSeqFile!=null){
		  writer=new BufferedWriter(new FileWriter(outSeqFile));			
		  for(SeqInfo seqInfo : seqObjList){		
			    seqIdentifier=seqInfo.seqIdentifier;
				seq=seqInfo.seq;
				seqQuality=seqInfo.seqQualityEncode;
				
				writer.write("@"+seqIdentifier);
				writer.newLine();
				writer.write(seq);
				writer.newLine();
				writer.write("+");
				writer.newLine();
				writer.write(seqQuality);
				writer.newLine();
				writer.flush();			
		  }
		  writer.close();
		  writer=null;
		}else{
		  System.err.println("null SeqInfo or outSeqFile!"); 
		}	
	}catch(Exception e){
		  System.out.println(e);
	}	 
  }
 
  public static void saveSeqObjAsFASTQ(SeqInfo seqInfo, String outSeqFile){
	 
	 try{
		String seqName="";
		String seq="";
		String seqQuality="";
		BufferedWriter writer=null;
		if(seqInfo!=null && outSeqFile!=null){
		    writer=new BufferedWriter(new FileWriter(outSeqFile));
		
		    seqName=seqInfo.seqName;
			seq=seqInfo.seq;
			seqQuality=seqInfo.seqQualityEncode;
			
			writer.write("@"+seqName+" "+seq.length());
			writer.newLine();
			writer.write(seq);
			writer.newLine();
			writer.write("+");
			writer.newLine();
			writer.write(seqQuality);
			writer.newLine();
			writer.flush();	
			
			writer.close();
			writer=null;
	   }else{
		  System.err.println("Empty SeqInfo or outSeqFile!"); 
	   }	
	}catch(Exception e){
		  System.out.println(e);
	}	 
  }  
  
  public static void saveSeqAsFASTA(String seq,String seqName,String outSeqFile){
	 
	 try{
	
		BufferedWriter writer=null;
		if(seq!=null && outSeqFile!=null){
		    writer=new BufferedWriter(new FileWriter(outSeqFile));		
            if(seqName==null) seqName="NA";
			writer.write(">"+seqName+" "+seq.length());
			writer.newLine();
			//writer.flush(); 
			writer.write(seq);
			writer.newLine();
			writer.flush();				
		
		    writer.close();
		    writer=null;
	   }else{
		  System.err.println("Empty SeqInfo or outSeqFile!"); 
	   }	
	}catch(Exception e){
		  System.out.println(e);
	}	 
  }
 
  public static void saveSeqObjAsFASTA(List<SeqInfo> seqObjList, BufferedWriter writer){
	 
	 try{
		String seqIdentifier="";
		String seq="";		
		if(seqObjList!=null && writer!=null){		
			
		  for(SeqInfo seqInfo : seqObjList){		
			seqIdentifier=seqInfo.seqIdentifier;
			seq=seqInfo.seq;
			writer.write(">"+seqIdentifier);
			writer.newLine();
			//writer.flush(); 
			writer.write(seq);
			writer.newLine();
			writer.flush();				
		  }
		}else{
		  System.err.println("null SeqInfo or outSeqFile!"); 
		}	
	}catch(Exception e){
		  System.out.println(e);
	}	 
  }
 
  public static void saveSeqObjAsFASTA(SeqInfo seqInfo, BufferedWriter writer){
	 
	 try{
		String seqName="";
		String seq="";	
		if(seqInfo!=null && writer!=null){		  
		
		    seqName=seqInfo.seqName;
			seq=seqInfo.seq;
			writer.write(">"+seqName+" "+seq.length());
			writer.newLine();
			//writer.flush(); 
			writer.write(seq);
			writer.newLine();
			writer.flush();		

	   }else{
		  System.err.println("Empty SeqInfo or outSeqFile!"); 
	   }	
	}catch(Exception e){
		  System.out.println(e);
	}	 
  }
 
  public static void saveSeqObjAsFASTQ(List<SeqInfo> seqObjList, BufferedWriter writer){
	 
	 try{
		String seqIdentifier="";
		String seq="";
		String seqQuality="";		
		if(seqObjList!=null && writer!=null){			
		  for(SeqInfo seqInfo : seqObjList){		
			    seqIdentifier=seqInfo.seqIdentifier;
				seq=seqInfo.seq;
				seqQuality=seqInfo.seqQualityEncode;
				
				writer.write("@"+seqIdentifier);
				writer.newLine();
				writer.write(seq);
				writer.newLine();
				writer.write("+");
				writer.newLine();
				writer.write(seqQuality);
				writer.newLine();
				writer.flush();			
		  }
		}else{
		  System.err.println("null SeqInfo or outSeqFile!"); 
		}	
	}catch(Exception e){
		  System.out.println(e);
	}	 
  }
 
  public static void saveSeqObjAsFASTQ(SeqInfo seqInfo, BufferedWriter writer){
	 
	 try{
		String seqName="";
		String seq="";
		String seqQuality="";
		
		if(seqInfo!=null && writer!=null){		
		    seqName=seqInfo.seqName;
			seq=seqInfo.seq;
			seqQuality=seqInfo.seqQualityEncode;
			
			writer.write("@"+seqName+" "+seq.length());
			writer.newLine();
			writer.write(seq);
			writer.newLine();
			writer.write("+");
			writer.newLine();
			writer.write(seqQuality);
			writer.newLine();
			writer.flush();	
	   }else{
		  System.err.println("Empty SeqInfo or outSeqFile!"); 
	   }	
	}catch(Exception e){
		  System.out.println(e);
	}	 
  }
 
  public static List<SeqInfo> combineSeqList( List<SeqInfo> seqList1, 
		 List<SeqInfo> seqList2){
	 
	 List<SeqInfo> seqObjList = new ArrayList<SeqInfo>();
	 for(SeqInfo seqInfo:seqList1){
		 seqObjList.add(seqInfo);
	 }
	 
	 for(SeqInfo seqInfo:seqList2){
		 seqObjList.add(seqInfo);
	 }
	 
	 return seqObjList;
  }

  public static List<String> combineStrList( List<String> strList1, 
		 List<String> strList2){
	 
	 List<String> strList = new ArrayList<String>();
	 for(String str:strList1){
		 if(!strList.contains(str))
		    strList.add(str);
	 }
	 
	 for(String str:strList2){
		 if(!strList.contains(str))
			    strList.add(str);
	 }
	 
	 return strList;
  }
 
  public static List<String> getGenomeChrLineSeq(String chrSeqFaFile){
		
		List<String> genomeChrLineSeqs;	
		genomeChrLineSeqs=loadChrSeqList(chrSeqFaFile);		
		
		return genomeChrLineSeqs;
		
  }
    
  public static List<String> loadChrSeqList(String chrSeqFaFile){
		
		List<String> chrSeqList=new ArrayList<String>();

		try{    
	        BufferedReader br;              
	        br = new BufferedReader(new FileReader(chrSeqFaFile));
		    String line;
			line = br.readLine();				
			while(true){           
		       if (line == null) break;
			   if(line.indexOf(">")>=0){			     
				 line=br.readLine();			
				 if (line == null) break;
				 while(line.indexOf(">")<0){
				   chrSeqList.add(line);
				   line=br.readLine();
				   if (line == null) break;
				 }				
			   }				 
			}           		   
			br.close();
			br=null;
		
		}catch(IOException e){
	        System.out.println(e);
	    }
				
	    return chrSeqList;
		
  }
    
  public static void setAlignSiteSeq(List<SeqAlignSite> chrAlignSites,List<String> genomeChrLineSeq,
		  int siteDownStreamLen,int siteUpStreamLen){
	  
	  for(ChrSite chrSite:chrAlignSites){
		  chrSite.seq=getChrSiteSeq(
				  genomeChrLineSeq,chrSite.chrStart,chrSite.chrEnd,
				  chrSite.strand,siteDownStreamLen,siteUpStreamLen
			  );
	  }
  }
 
  public static void setChrSiteSeq(List<ChrSite> chrSites,List<String> genomeChrLineSeq,
		  int siteDownStreamLen,int siteUpStreamLen){
	  
	  if(chrSites==null) return;
	  
	  for(ChrSite chrSite:chrSites){
		  chrSite.seq=getChrSiteSeq(
			  genomeChrLineSeq,chrSite.chrStart,chrSite.chrEnd,
			  chrSite.strand,siteDownStreamLen,siteUpStreamLen
		  );
	  }
  }
 
  public static void setChrSiteSeq(List<ChrSite> chrSites,List<String> genomeChrLineSeq,
		  int siteDownStreamLen,int siteUpStreamLen, boolean doRCForMinusStrand){
	  
	  if(chrSites==null) return;
	  if(doRCForMinusStrand) {
	    for(ChrSite chrSite:chrSites){
		  chrSite.seq=getChrSiteSeq(
				  genomeChrLineSeq,chrSite.chrStart,chrSite.chrEnd,
				  chrSite.strand,siteDownStreamLen,siteUpStreamLen
		  );
		  if(chrSite.strand.equals("-")) chrSite.seq=SeqOperation.reverseComplement(chrSite.seq);
	    }
	  }else {
		for(ChrSite chrSite:chrSites){
		    chrSite.seq=getChrSiteSeq(
				  genomeChrLineSeq,chrSite.chrStart,chrSite.chrEnd,
				  chrSite.strand,siteDownStreamLen,siteUpStreamLen
			);
	    }
	  }
  }
  
  public static String getChrSiteSeq(List<String> genomeChrLineSeq,int siteStart,int siteEnd,
		  String siteStrand,int siteDownStreamLen,int siteUpStreamLen){

		String siteDownStreamSeq = null;
		String siteUpStreamSeq = null;
		String siteAroundSeq = null;
		
		int coveredDownFaRow;
		int coveredUpFaRow;
		int seqRowStart;
		int seqRowEnd;
		int seqStart;
		//int seqEnd;
		int tmpLength;
		
		int seqSizePerLine=genomeChrLineSeq.get(0).length();		
	
		siteDownStreamLen=(siteEnd-siteStart)+siteDownStreamLen;
		coveredDownFaRow=siteDownStreamLen/seqSizePerLine+1;
		coveredUpFaRow=siteUpStreamLen/seqSizePerLine+1;	
		
		//System.out.println("GC Content for "+chrName);	  	
		int start=siteStart;
		if(siteStrand.equals("+")) 
		   start=siteStart;
		else if(siteStrand.equals("-")) 
		   start=siteEnd;
				
		double eventFaRowCoord;
		int eventFaRowCoord_iPart;
		////double eventFaRowCoord_fPart;
		int eventFaRowCoord_fPart_pos;
		
		eventFaRowCoord=(1.0d*start)/(seqSizePerLine*1.0d);
		eventFaRowCoord_iPart=(int) Math.floor(eventFaRowCoord);
		////eventFaRowCoord_fPart=eventFaRowCoord-Math.floor(eventFaRowCoord);
		////eventFaRowCoord_fPart_pos=(int) (seqSizePerLine*eventFaRowCoord_fPart);
		eventFaRowCoord_fPart_pos=start%seqSizePerLine;
		
		if(siteStrand.equals("+")){			  
		  siteUpStreamSeq="";
		  seqRowStart=eventFaRowCoord_iPart-coveredUpFaRow;
		  if(seqRowStart<0) seqRowStart=0;
		  for(int n=seqRowStart;n<eventFaRowCoord_iPart;n++){
			siteUpStreamSeq=siteUpStreamSeq+genomeChrLineSeq.get(n);
		  }
		  siteUpStreamSeq=siteUpStreamSeq
				+genomeChrLineSeq.get(eventFaRowCoord_iPart).substring(0,eventFaRowCoord_fPart_pos);
		  seqStart=siteUpStreamSeq.length()-siteUpStreamLen-1;
		  if(seqStart<0) seqStart=0;
		  siteUpStreamSeq=siteUpStreamSeq.substring(seqStart,siteUpStreamSeq.length());
		  
		  tmpLength=genomeChrLineSeq.get(eventFaRowCoord_iPart).length();
		  siteDownStreamSeq=genomeChrLineSeq.get(eventFaRowCoord_iPart).substring(
			   eventFaRowCoord_fPart_pos,Math.min(tmpLength,seqSizePerLine)
		  );
		  seqRowEnd=eventFaRowCoord_iPart+1+coveredDownFaRow;
		  if(seqRowEnd>genomeChrLineSeq.size()) seqRowEnd=genomeChrLineSeq.size();
		  for(int n=eventFaRowCoord_iPart+1;n<seqRowEnd;n++){
			siteDownStreamSeq=siteDownStreamSeq+genomeChrLineSeq.get(n);
		  }
			
		  if(siteDownStreamSeq.length()>siteDownStreamLen) {
			siteDownStreamSeq=siteDownStreamSeq.substring(0,siteDownStreamLen);
		  }else {
			siteDownStreamSeq=siteDownStreamSeq.substring(0,siteDownStreamSeq.length());	
		  }			
		  siteAroundSeq=siteUpStreamSeq+siteDownStreamSeq;		 
		}else if(siteStrand.equals("-")){				  
		  siteDownStreamSeq="";
		  seqRowStart=eventFaRowCoord_iPart-coveredDownFaRow;
		  if(seqRowStart<0) seqRowStart=0;
		  for(int n=seqRowStart;n<eventFaRowCoord_iPart;n++){
			siteDownStreamSeq=siteDownStreamSeq+genomeChrLineSeq.get(n);
		  }
		  siteDownStreamSeq=siteDownStreamSeq
				+genomeChrLineSeq.get(eventFaRowCoord_iPart).substring(0,eventFaRowCoord_fPart_pos);
		  seqStart=siteDownStreamSeq.length()-siteDownStreamLen-1;
		  if(seqStart<0) seqStart=0;
		  siteDownStreamSeq=siteDownStreamSeq.substring(seqStart,siteDownStreamSeq.length());
					
		  siteUpStreamSeq="";
		  tmpLength=genomeChrLineSeq.get(eventFaRowCoord_iPart).length();
		  if(tmpLength>=seqSizePerLine) {/////
			 siteUpStreamSeq=
				genomeChrLineSeq.get(eventFaRowCoord_iPart).substring(eventFaRowCoord_fPart_pos,seqSizePerLine);
		  }else{
			siteUpStreamSeq=
				genomeChrLineSeq.get(eventFaRowCoord_iPart).substring(eventFaRowCoord_fPart_pos,tmpLength);
		  }/////
			
		  seqRowEnd=eventFaRowCoord_iPart+1+coveredUpFaRow;
		  if(seqRowEnd>genomeChrLineSeq.size()) seqRowEnd=genomeChrLineSeq.size();
		  for(int n=eventFaRowCoord_iPart+1;n<seqRowEnd;n++){
			siteUpStreamSeq=siteUpStreamSeq+genomeChrLineSeq.get(n);
		  }
		  if(siteUpStreamSeq.length()>siteUpStreamLen) {
			siteUpStreamSeq=siteUpStreamSeq.substring(0,siteUpStreamLen);
		  }else {
			siteUpStreamSeq=siteUpStreamSeq.substring(0,siteUpStreamSeq.length());	
		  }					    
		  siteAroundSeq=siteDownStreamSeq+siteUpStreamSeq;			
		}
		
		return siteAroundSeq;
		
  }  

}