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
package org.geatools.call;

import java.io.File;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Calendar;
import java.util.List;
import java.util.Map;

import org.geatools.GEAT;
import org.geatools.data.structure.SeqMutInfo;
import org.geatools.data.structure.SeqQual;
import org.geatools.data.structure.SeqRocket;
import org.geatools.data.structure.SeqRocketPair;
import org.geatools.operation.FileOperate;
import org.geatools.seqprocess.SeqMut;
import org.geatools.seqprocess.SeqOperation;
import org.geatools.seqprocess.SeqQCFilter;
import org.geatools.seqprocess.SeqRocketConsole;

public class CallSeqMut extends GEAT{
	static String homeDir=null;	
	static String workingDir=null;
	static String tmpDir=null;
	
	//Alignment Quality Standard(AQS) filter criteria
	static int minAlignLen=100;
	static int maxMismatch=9999;
	static int maxGapNum=9999;
	static int maxQStart=50;
	static int minQStart=1;
	static int maxQEnd=500;
	static int minQEnd=200;
	static int maxSStart=50;
	static int minSStart=1;
	static int maxSEnd=500;
	static int minSEnd=200;
	
	//Neighbor Quality Standard(NQS) filter criteria
	static int NQS_minS=30;
	static int NQS_NN=5;
	static int NQS_minNS=25;
	
	//INDEL filter criteria
	static int maxIndel=1;
	
	//Local Mutation Standard(LMS) filter criteria
	static int localBaseNum=5;
	static int localMaxMutNum=3;
	
	//base mut rate filter criteria
	static float snpRate=0.2f;
	static float baseErrRate=0.01f;
	
	//for mut of particular bases, for example C->T, G->A 
	static String [][]interestMut;
	
	static boolean skipAQS=true;
	static boolean skipNQS=true;
	static boolean skipIndel=true;
	static boolean skipLMS=true;
	static boolean skipBaseMutCheck=true;
	static boolean doSNPCheck=false;
	static boolean doBaseErrCheck=false;
	
	static int avgRegion_start=1;
	static int avgRegion_end=200;
	
	public static void setHomeDir(String dir){
		homeDir=dir;
	}
	public static void setWorkingDir(String dir){
		workingDir=dir;
	}
	public static void setTmpDir(String dir){
		tmpDir=dir;
	}
	public static void delTmpDir(String dir){
		FileOperate.delFolder(dir);
	}
	
	public static void doWork(String[] args) {
		  
		 String querySeqFile = null;
		 String querySeqFile2=null;
		 String querySeqListFile = null;
		 String querySeqListFile2 = null;
		 String refSeqFile = null;
		 List<String> seqQCOptList;
		 List<String>querySeqList; 
		 List<String>querySeqList2; 
		 List<SeqQual> seqQualList=null;
		 List<SeqQual> seqQualList2=null;		
		 
		 if(homeDir==null) homeDir=GEAT.getClassPath();			   
		 if(tmpDir==null){
			 String timeStamp=new SimpleDateFormat("yyyyMMdd_HHmmss").format(Calendar.getInstance().getTime());
			 tmpDir=homeDir+"/tmp/"+timeStamp;
		 } 
		 if(workingDir==null) workingDir=homeDir+"/working";
		 File dir=new File(tmpDir);
		 if(!dir.exists()) FileOperate.newFolder(tmpDir);
		 dir=null;
		 dir=new File(workingDir);
		 if(!dir.exists()) FileOperate.newFolder(workingDir);
		 dir=null;
			
		 Map<String, List<String>> params=getCommandLineParams(args);
		 if(params.get("task")!=null){			
			if(params.get("task").size()>0){
			  taskName=params.get("task").get(0).trim();
			  if(taskName!=null){
			    if(!taskName.equalsIgnoreCase("SeqMut")){			      
			      System.err.println("Error: '-task' is invalid");
				  return;	
			    }
			  }
			}
		 }else{
		   System.err.println("Error: '-task' is invalid");
		   return;		
		 }
		 
		 if(params.get("fastq")!=null){
			 isFastqOK=false;
			 if(params.get("fastq").size()>0){
		       fastq=params.get("fastq").get(0).trim();
		       if(fastq!=null){
		    	 if(SeqOperation.isFASTQSeq(fastq)){
		    	   isFastqOK=true;
		    	   //doSeqQCFilter=true;
		    	   seqType=SeqOperation.SEQTYPE_SINGLEEND;
		    	 }else{
			       System.err.println("The -fastq file doesn't exist or isn't a fastq file:(");
			       return;
			     }		    	
		       }
		     }else{
		       System.err.println("You didn't provide fastq file :(");
		       return;
		     }
		  }else if(params.get("fasta")!=null){
			  isFastaOK=false;
			  if(params.get("fasta").size()>0){
				 fasta=params.get("fasta").get(0).trim();
				 if(fasta!=null){
					 if(SeqOperation.isFASTASeq(fasta)){
				    	isFastaOK=true;
				    	doSeqQCFilter=false;
				    	seqType=SeqOperation.SEQTYPE_SINGLEEND;				    	
				     }else{
					    System.err.println("The -fasta file doesn't exist or isn't a fasta file:(");
					    return;
					 }				 
				  }
			   }else{
				  System.err.println("You didn't provide fasta file :(");
				  return;
			   }
		   }
		   
		   if(params.get("fastq2")!=null){
			   isFastq2OK=false;
			   if(params.get("fastq2").size()>0){
			       fastq2=params.get("fastq2").get(0).trim();
			       if(fastq2!=null){
			    	 if(SeqOperation.isFASTQSeq(fastq2)){
			    	   isFastq2OK=true;
			    	   //doSeqQCFilter=true;
			    	   seqType=SeqOperation.SEQTYPE_PAIREND;
			    	 }else{
				       System.err.println("The -fastq2 file doesn't exist or isn't a fastq file :(");
				       return;
				     }			    
			       }
			   }else{
			       System.err.println("You didn't provide fastq2 file :(");
			       return;
			   }
		   }else if(params.get("fasta2")!=null){
			   isFasta2OK=false;
			   if(params.get("fasta2").size()>0){
				   fasta2=params.get("fasta2").get(0).trim();
				   if(fasta2!=null){
					  if(SeqOperation.isFASTASeq(fasta2)){
					    isFasta2OK=true;
					    doSeqQCFilter=false;
					    seqType=SeqOperation.SEQTYPE_PAIREND;					   
					  }else{
						System.err.println("The -fasta2 file doesn't exist or isn't a fasta file :(");
						return;
					  }				
				   }
				}else{
				   System.err.println("You didn't provide fasta2 file :(");
				   return;
				}
		   }
		  
		   //####### set querySeq #######
		   querySeqList=new ArrayList<String>();
		   if(params.get("querySeq")!=null){		    
		      if(params.get("querySeq").size()>0){
		    	 
				 querySeqFile=params.get("querySeq").get(0).trim();					
				 if(SeqOperation.isFASTASeq(querySeqFile) 
						 || SeqOperation.isFASTQSeq(querySeqFile)){					
					
					querySeqList.add(querySeqFile);
					seqType=SeqOperation.SEQTYPE_SINGLEEND;
					if(SeqOperation.isFASTQSeq(querySeqFile)) doSeqQCFilter=true;
					 
				 }else{
					System.err.println("Error: '-querySeq' file doesn't exist or isn't a fasta/fastq file.");
					querySeqFile=null;
					return;
				 }					 
			  }else{
				 System.err.println("Illegal '-querySeq' parameter usage :(");				
				 return;
			  }
		   }
		   
		   querySeqList2=new ArrayList<String>();
		   if(params.get("querySeq2")!=null){		    	
		      if(params.get("querySeq2").size()>0){
				 querySeqFile2=params.get("querySeq2").get(0).trim();					
				 if(SeqOperation.isFASTASeq(querySeqFile2) 
						 || SeqOperation.isFASTQSeq(querySeqFile2)){
					
					querySeqList2.add(querySeqFile2);
					seqType=SeqOperation.SEQTYPE_PAIREND;	
					if(SeqOperation.isFASTQSeq(querySeqFile2)) doSeqQCFilter=true;
				 
				 }else{ 
					System.err.println("Error: '-querySeq2' file doesn't exist or isn't a fasta/fastq file.");
					querySeqFile2=null;					
					return;
				 }	
			  }else{
				 System.err.println("Illegal '-querySeq2' parameter usage :(");		
				 return;
			  }
		   }
		   
		   //####### set querySeq #######
		   if(params.get("querySeqList")!=null){		    	
		     if(params.get("querySeqList").size()>0){
				querySeqListFile=params.get("querySeqList").get(0).trim();
				querySeqList=FileOperate.getRowListFromFile(querySeqListFile);
				if(SeqOperation.isFASTASeq(querySeqList) 
						 || SeqOperation.isFASTQSeq(querySeqList)){
					
					seqType=SeqOperation.SEQTYPE_SINGLEEND;
					if(SeqOperation.isFASTQSeq(querySeqList)) doSeqQCFilter=true;
					
				}else{	
					System.err.println("Error: '-querySeqList' file doesn't exist or doesn't contain a fasta/fastq file.");
					querySeqListFile=null;
					querySeqList=new ArrayList<String>();
					return;
				}					 
			 }else{
				System.err.println("Illegal '-querySeqList' parameter usage :(");			
				return;
			 }
		   }
		   
		   if(params.get("querySeqList2")!=null){		    	
		      if(params.get("querySeqList2").size()>0){
				 querySeqListFile2=params.get("querySeqList2").get(0).trim();					
				 querySeqList2=FileOperate.getRowListFromFile(querySeqListFile2);
				 if(SeqOperation.isFASTASeq(querySeqList2) 
						 || SeqOperation.isFASTQSeq(querySeqList2)){
					
					seqType=SeqOperation.SEQTYPE_PAIREND;
					if(SeqOperation.isFASTQSeq(querySeqList2)) doSeqQCFilter=true;
				 
				 }else{						
					System.err.println("Error: '-querySeqList2' file doesn't exist or doesn't contain a fasta/fastq file.");
					querySeqListFile2=null;	
					querySeqList2=new ArrayList<String>();
				    return;				    
				 }				
			  }else{
				 System.err.println("Illegal '-querySeqList2' parameter usage :(");			
				 return;
			  }
		   }
		     
		   //####### set refSeq #######
		   if(params.get("refSeq")!=null){		    
		      if(params.get("refSeq").size()>0){
				 refSeqFile=params.get("refSeq").get(0).trim();					
				 if(!SeqOperation.isFASTASeq(refSeqFile)){
					System.err.println("Error: '-refSeq' file doesn't exist or isn't a fasta file.");						
					return;
				 }					
			  }else{
				 System.err.println("Illegal '-refSeq' parameter usage :(");				
				 return;
			  }
		   }
		   //check -querySeq, -refSeq
		   if((querySeqFile==null && !isFastqOK && !isFastaOK) && (querySeqListFile==null)){
			  System.err.println("Error: missed/invalid '-querySeq or -fasta or -fastq' parameter!!!");
		
			  return; 
		   }
		   if(refSeqFile==null){
			  System.err.println("Error: missed/invalid '-refSeq' parameter!!!");
			
			  return; 
		   }
		   
		   //####### set output #######
		   boolean  doOutput=false;
		   if(params.get("outName")!=null){				
				if(params.get("outName").size()>0){
					 outName=params.get("outName").get(0).trim();					
				}else{
					 System.err.println("Illegal '-outName' parameter usage :(");
					 return;
				}
		   }
			
		   if(params.get("outTag")!=null){		  	
				if(params.get("outTag").size()>0){
					 outTag=params.get("outTag").get(0).trim();						 
				}else{
					 System.err.println("Illegal '-outTag' parameter usage :(");
					 return;
				}
		   }
			
		   if(params.get("outDir")!=null){			
			   if(params.get("outDir").size()>0){
				 outDir=params.get("outDir").get(0).trim();
				 File f=new File(outDir);
			     if (f.exists()){
				   doOutput=true;
				   f=null;
				 }else if (outDir!=null){
				   FileOperate.newFolder(outDir);
				   doOutput=true;
				 }			  
			   }else{
				 System.err.println("Illegal '-outDir' parameter usage :(");		
				 return;
				}
		   }		   
		   if(!doOutput){			   
				 outDir=homeDir+"/working";
				 dir=new File(outDir);
				 if(!dir.exists()) FileOperate.newFolder(outDir);	
		         dir=null;
				 doOutput=true;   			
		   }
		   
		   //####### set seq info for experiments in library#######
		   if(params.get("expSeqInfo")!=null){
			  isLibExpSeqInfoOK=false;
			  if(params.get("expSeqInfo").size()>0){
			     expSeqInfo=params.get("expSeqInfo").get(0).trim();
			     if(expSeqInfo!=null){
			       File f = new File(expSeqInfo);
			       if(f.exists()){
			    	 isLibExpSeqInfoOK=true;
			    	 seqType=SeqOperation.SEQTYPE_SINGLEEND;
			       }else{
				     System.err.println(" '-expSeqInfo' file you provided doesn't exist :(");
				     return;
				   }
			       f=null;
			     }
			  }else{
			     System.err.println("You didn't provide '-expSeqInfo' file :(");
			     return;
			  }
		   }
		   
		   if(params.get("expSeqInfo2")!=null){
			  isLibExpSeqInfoOK=false;
			  if(params.get("expSeqInfo2").size()>0){
			     expSeqInfo2=params.get("expSeqInfo2").get(0).trim();
			     if(expSeqInfo2!=null){
			       File f = new File(expSeqInfo2);
			       if(f.exists()){
			    	 isLibExpSeqInfoOK=true;
			    	 seqType=SeqOperation.SEQTYPE_PAIREND;
			       }else{
				     System.err.println(" '-expSeqInfo2' file you provided doesn't exist :(");
				     return;
				   }
			       f=null;
			     }
			  }else{
			     System.err.println("You didn't provide '-expSeqInfo2' file :(");
			     return;
			  }
		   }
		   
		   seqQCOptList=new ArrayList<String>();		
		   if(params.get("filter_seqQC")!=null){		
			   doSeqQCFilter=true;
			   String subParams;
			   String [] itemSplited;
			   if(params.get("filter_seqQC").size()>0){
				 subParams=params.get("filter_seqQC").get(0).trim();
				 if(subParams!=null && subParams.equalsIgnoreCase("skip"))			    		  
				    doSeqQCFilter=false;
			     else{	
				    //List<String> seqQCOptList=new ArrayList<String>();
				    for(int i=0;i<params.get("filter_seqQC").size();i++){				
					  subParams=params.get("filter_seqQC").get(i).trim();
					  itemSplited=subParams.split("=");
					  if(itemSplited.length>1){
						if(itemSplited[0].trim().equalsIgnoreCase("min_qual_mean")){
							seqQCOptList.add("-min_qual_mean");
							seqQCOptList.add(itemSplited[1].trim());
						}else if(itemSplited[0].trim().equalsIgnoreCase("min_len")){
							seqQCOptList.add("-min_len");
							seqQCOptList.add(itemSplited[1].trim());
						}else if(itemSplited[0].trim().equalsIgnoreCase("de_exact_dup")){
							seqQCOptList.add("-de_exact_dup");
							seqQCOptList.add(itemSplited[1].trim());
						}					
					  }
				    }
			      }
			   }else{
				  System.out.println("Warning: empty '-filter_seqQC' parameter,the system uses default '-filter_seqQC' parameter.");
				  doSeqQCFilter=true;
			   }				
		   }
		   
		   if(params.get("filter_AQS")!=null){		
			   skipAQS=false;
			   String subParams;
			   String [] itemSplited;
			   if(params.get("filter_AQS").size()>0){
				 subParams=params.get("filter_AQS").get(0).trim();
				 if(subParams!=null && subParams.equalsIgnoreCase("skip")) skipAQS=true;						 
			     else{				
				   for(int i=0;i<params.get("filter_AQS").size();i++){				
					 subParams=params.get("filter_AQS").get(i).trim();
					 itemSplited=subParams.split("=");
					 if(itemSplited.length>1){
						try{
						  if(itemSplited[0].trim().equalsIgnoreCase("minAlignLen")){
							minAlignLen=Integer.parseInt(itemSplited[1].trim());
						  }else if(itemSplited[0].trim().equalsIgnoreCase("maxMismatch")){
							maxMismatch=Integer.parseInt(itemSplited[1].trim());
						  }else if(itemSplited[0].trim().equalsIgnoreCase("maxGapNum")){
							maxGapNum=Integer.parseInt(itemSplited[1].trim());
						  }
						}catch(Exception ex){
						  System.out.println("Warning: The system will use default '-filter_AQS' parameter.");
						  skipAQS=false;	
						}
					 }
				   }
			     }
			   }else{
				 System.out.println("Warning: The system will use default '-filter_AQS' parameter.");
				 skipAQS=false;	
			   }				
		   }
		   
		   if(params.get("filter_NQS")!=null){		
			   skipNQS=false;	
			   String subParams;
			   String [] itemSplited;
			   if(params.get("filter_NQS").size()>0){
				 subParams=params.get("filter_NQS").get(0).trim();
				 if(subParams!=null && subParams.equalsIgnoreCase("skip")) skipNQS=true;						 
			     else{						
				    for(int i=0;i<params.get("filter_NQS").size();i++){				
					  subParams=params.get("filter_NQS").get(i).trim();
					  itemSplited=subParams.split("=");
					  if(itemSplited.length>1){
						try{ 
						  if(itemSplited[0].trim().equalsIgnoreCase("minScore")){
							NQS_minS=Integer.parseInt(itemSplited[1].trim());
						  }else if(itemSplited[0].trim().equalsIgnoreCase("minNScore")){
							NQS_minNS=Integer.parseInt(itemSplited[1].trim());
						  }else if(itemSplited[0].trim().equalsIgnoreCase("NNum")){
							NQS_NN=Integer.parseInt(itemSplited[1].trim());
						  }
					    }catch(Exception ex){
						  System.out.println("Warning: The system will use default '-filter_NQS' parameter.");
						  skipNQS=false;	
					    }	
					  }
				    }
			     }
			   }else{
				  System.out.println("Warning: The system will use default '-filter_NQS' parameter.");
				  skipNQS=false;	
			   }				
		   }  
		   
		   if(params.get("filter_INDEL")!=null){		
			   skipIndel=false;	
			   String subParams;
			   String [] itemSplited;
			   if(params.get("filter_INDEL").size()>0){
				 subParams=params.get("filter_INDEL").get(0).trim();
				 if(subParams!=null && subParams.equalsIgnoreCase("skip")) skipIndel=true;						 
			     else{				
				    for(int i=0;i<params.get("filter_INDEL").size();i++){				
					  subParams=params.get("filter_INDEL").get(i).trim();
					  itemSplited=subParams.split("=");
					  if(itemSplited.length>1){
						try{
						  if(itemSplited[0].trim().equalsIgnoreCase("maxIndel")){
							maxIndel=Integer.parseInt(itemSplited[1].trim());
						  }
						}catch(Exception ex){
						  System.out.println("Warning: The system will use default '-filter_INDEL' parameter.");
						  skipIndel=false;	
						}
					  }
				    }
			     }
			   }else{
				  System.out.println("Warning: The system will use default '-filter_INDEL' parameter.");
				  skipIndel=false;	
			   }				
		   }  
		   
		   if(params.get("filter_LMS")!=null){		
			   skipLMS=false;	
			   String subParams;
			   String [] itemSplited;
			   if(params.get("filter_LMS").size()>0){
				 subParams=params.get("filter_LMS").get(0).trim();
				 if(subParams!=null && subParams.equalsIgnoreCase("skip")) skipLMS=true;						 
			     else{
				   for(int i=0;i<params.get("filter_LMS").size();i++){				
					  subParams=params.get("filter_LMS").get(i).trim();
					  itemSplited=subParams.split("=");
					  if(itemSplited.length>1){
						try{
						  if(itemSplited[0].trim().equalsIgnoreCase("baseNum")){
							 localBaseNum=Integer.parseInt(itemSplited[1].trim());
						  }else if(itemSplited[0].trim().equalsIgnoreCase("maxMutNum")){
							 localMaxMutNum=Integer.parseInt(itemSplited[1].trim());
						  }
						}catch(Exception ex){
						  System.out.println("Warning: The system will use default '-filter_LMS' parameter.");
						  skipLMS=false;	
						}
					  }
				   }
			     }
			   }else{
				  System.out.println("Warning: The system will use default '-filter_LMS' parameter.");
				  skipLMS=false;	
			   }				
		   }
		   
		   if(params.get("filter_BMF")!=null){	
			   skipBaseMutCheck=false;
			   doSNPCheck=false;
			   doBaseErrCheck=false;
			   String subParams;
			   String [] itemSplited;
			   if(params.get("filter_BMF").size()>0){
				 subParams=params.get("filter_BMF").get(0).trim();
				 if(subParams!=null && subParams.equalsIgnoreCase("skip")){
				    skipBaseMutCheck=true;
					doSNPCheck=false;
					doBaseErrCheck=false;
				 }else{
				   for(int i=0;i<params.get("filter_BMF").size();i++){				
					  subParams=params.get("filter_BMF").get(i).trim();
					  itemSplited=subParams.split("=");
					  if(itemSplited.length>1){
						try{
						  if(itemSplited[0].trim().equalsIgnoreCase("snpRate")){
							 snpRate=Float.parseFloat(itemSplited[1].trim());
							 doSNPCheck=true;
						  }else if(itemSplited[0].trim().equalsIgnoreCase("errRate")){
							 baseErrRate=Float.parseFloat(itemSplited[1].trim());
							 doBaseErrCheck=true;
						  }
						}catch(Exception ex){
						  System.out.println("Warning: The system will use default '-filter_BMF' parameter.");
						}
					  }
				   }
			     }
			   }else{
				  System.out.println("Warning: The system will use default '-filter_BMF' parameter.");
				  skipBaseMutCheck=false;
				  doSNPCheck=true;
				  doBaseErrCheck=true;
				  snpRate=0.2f;
				  baseErrRate=0.01f;
			   }				
		   }
		   
		   if(params.get("tarMut")!=null){		
			 if(params.get("tarMut").size()>0){		
				String subParams;
				String [] itemSplited;
				List<String> refBases=new ArrayList<String>();
				List<String> altBases=new ArrayList<String>();
				for(int i=0;i<params.get("tarMut").size();i++){
				   subParams=params.get("tarMut").get(i).trim();
				   itemSplited=subParams.split("2");			
				   if(itemSplited.length>1 
						   && itemSplited[0].length()==1
						   && itemSplited[1].length()==1){				 
				     refBases.add(itemSplited[0]);
				     altBases.add(itemSplited[1]);
				   }
				}
				
				if(refBases.size()>0){
				   interestMut=new String[refBases.size()][2];
				   for(int i=0;i<interestMut.length;i++){
					   interestMut[i][0]=refBases.get(i);
					   interestMut[i][1]=altBases.get(i);
				   }
				   //refBases=null;
				   //altBases=null;
				}else{
				   System.out.println("Warning: '-tarMut' parameter is illegal. The system will skip it.");			
				}
			 }else{
				System.out.println("Warning: '-tarMut' parameter is empty. The system will skip it.");			
			 }
				  
		   }
		   
		   if(params.get("avgRegion")!=null){
			   String subParams;
			   String [] itemSplited;
			   if(params.get("avgRegion").size()>0){		
				   for(int i=0;i<params.get("avgRegion").size();i++){				
					  subParams=params.get("avgRegion").get(i).trim();
					  itemSplited=subParams.split("=");
					  if(itemSplited.length>1){
						try{
						  if(itemSplited[0].trim().equalsIgnoreCase("start")){
							  avgRegion_start=Integer.parseInt(itemSplited[1].trim());
						  }else if(itemSplited[0].trim().equalsIgnoreCase("end")){
							  avgRegion_end=Integer.parseInt(itemSplited[1].trim());
						  }
						  
						  if(avgRegion_end<avgRegion_start){
							System.err.println("Error: 'end' must be larger than 'start' for '-avgRegion' parameter.");
							return;
						  }
						}catch(Exception ex){
						  System.out.println("Warning: The system will use default '-avgRegion' parameter.");				
						}
					  }
				   }			     
			   }else{
				  System.out.println("Warning: The system will use default '-avgRegion' parameter.");	
			   }				
		   }
		   
		   //####### split seq from given  fasta or fastq file########
		   String splitedSeqOut=null;
		   if(params.get("splitSeq")!=null){
			 doSplitSeq=true;				
			 String subParams;
			 String [] itemSplited;
			 if(params.get("splitSeq").size()>0){
			   subParams=params.get("splitSeq").get(0).trim();
			   if(subParams!=null && subParams.equalsIgnoreCase("skip")) doSplitSeq=false;						 
			   else{
				  for(int i=0;i<params.get("splitSeq").size();i++){				
					 subParams=params.get("splitSeq").get(i).trim();
					 itemSplited=subParams.split("=");
					 if(itemSplited.length>1){
						if(itemSplited[0].equalsIgnoreCase("step")){
						   if(isInteger(itemSplited[1])){ 
							 splitStep=Integer.parseInt(itemSplited[1]); 
						   }else{
							 System.out.println("Warning: '-splitSeq step' parameter is illeagl.");
							 System.out.println("Warning: the system will use default '-splitSeq step' parameters.");
							 doSplitSeq=true;
						   }
						}else if(itemSplited[0].equalsIgnoreCase("out")){
						   splitedSeqOut=itemSplited[1];
						}
					 }
				  }
			   }
			 }else{
			   System.out.println("Warning: the system will use default '-splitSeq' parameters.");
			   doSplitSeq=true;
			 }
		   }
		   
		   //for forward seq
		   //-fastq or -fasta as forward raw seq for recognizing valid sample sequences or demultiplexing samples
		   //It must be provided together with -expSeqInfo which contains info for recognizing valid sample sequences or demultiplexing samples
		   //if if both fastq and fasta are provided, we only consider fastq. 
		   String inSeqFile="";
		   if(isFastqOK && fastq!=null){
			  inSeqFile=fastq;			 
		   }else if(isFastaOK && fasta!=null){	
	    	  inSeqFile=fasta;	    	  
	       }
		   
		   // for reverse seq
		   //-fastq2 or -fasta2 as reverse raw seq for recognizing valid sample sequence or demultiplexing samples
		   //It must be provided together with -expSeqInfo2 which contains info for recognizing valid sample sequence or demultiplexing samples
		   //if if both fastq2 and fast2a are provided, we only consider fastq2. 
		   String inSeqFile2="";
		   if(isFastq2OK && fastq2!=null){
			  inSeqFile2=fastq2;
		   }else if(isFasta2OK && fasta2!=null){	
	    	  inSeqFile2=fasta2;
	       }		   
		
		   if(isFastqOK || isFastaOK || isFastq2OK || isFasta2OK){
			   if(!isLibExpSeqInfoOK){
				  System.err.println("Error, missing or invalid -expSeqInfo parameter!");
				  return; 
			   }
		   }
		   
		   //=======================start to split seq===========================
		   if(doSplitSeq &&(isFastqOK || isFastaOK)){
				splitedSeqFiles=null;
				//Check seq format, and then split seq into multiple subfiles................
				System.out.println(".........Spliting forward seq..........");
				splitedSeqOut=outDir+"/split_forward";
				FileOperate.newFolder(splitedSeqOut);
	      	    splitedSeqFiles=SeqOperation.splitSeqFile(inSeqFile,splitStep,
	      	    		splitedSeqOut);
		 				
				if(splitedSeqFiles==null || splitedSeqFiles.size()==0){
					System.err.println("Sorry, not working for sequences split!");
				    return;
				}
				
				if(seqType.equalsIgnoreCase(SeqOperation.SEQTYPE_PAIREND)){				
					splitedSeqFiles2=null;	
					//Check seq format, and then split seq into multiple subfiles................
					System.out.println(".........Spliting pair-end reverse seq.........");
					splitedSeqOut=outDir+"/split_reverse";
					FileOperate.newFolder(splitedSeqOut);
			        splitedSeqFiles2=SeqOperation.splitSeqFile(inSeqFile2,splitStep,
			        		splitedSeqOut);			
			    	
			    	if(splitedSeqFiles==null || splitedSeqFiles.size()==0){
					   System.err.println("Sorry, not working for sequences split!");
				       return;
				    }
				}				
		   }	
		   
		   //===================Recognizing sequence for each experiment=================
		   //Demultiplexing samples from given fastq librrary
		   if((isFastqOK || isFastaOK) && isLibExpSeqInfoOK){
	    	   System.out.println("..........Recognizing sequence for each experiment..........");
	    	   isSeqRocketOK=false;
	    	   seqRC=new SeqRocketConsole();    	   
	    	   seqRC.setHomeDir(homeDir);
	    	   seqRC.setDataDir(dataDir);
	    	   seqRC.setTmpDir(tmpDir);
	    	   
	    	   if(seqType.equalsIgnoreCase(SeqOperation.SEQTYPE_SINGLEEND)){
	    		 if(doSplitSeq)
	    		   seqRC.splitLaunchSingleEnd(splitedSeqFiles,expSeqInfo,outDir);
	    		 else
	    		   seqRC.launchSingleEnd(inSeqFile,expSeqInfo,outDir);
	    		 
			       isSeqRocketOK=seqRC.isSeqRocketsOK();
				 
			     //set querySeqList by demultiplexing samples,querySeqList2 for reverse
				 if(isSeqRocketOK){ 				
					  querySeqList=new ArrayList<String>();
					  for(SeqRocket seqRocket: seqRC.getSeqRockets()){							  
						querySeqList.add(seqRocket.recognizedSeqFile);
					  }
				  }

	    	   }else if(seqType.equalsIgnoreCase(SeqOperation.SEQTYPE_PAIREND)){
	    		   if(doSplitSeq)
	    		     seqRC.splitLaunchPairEnd(splitedSeqFiles,splitedSeqFiles2,
			        		expSeqInfo,expSeqInfo2,outDir);	
	    		   else
	    			 seqRC.launchPairEnd(inSeqFile,inSeqFile2,expSeqInfo,expSeqInfo2,
	    					 outDir);
	    		   
		    	   isSeqPairRocketOK=seqRC.isSeqRocketsOK();
		    	   
		    	 //set querySeqList by demultiplexing samples,querySeqList2 for reverse
				   if(isSeqPairRocketOK){
					  querySeqList=new ArrayList<String>();
					  querySeqList2=new ArrayList<String>();
					  for(SeqRocketPair seqPairRocket:seqRC.getSeqPairRockets()){	
						 querySeqList.add(seqPairRocket.forward.recognizedSeqFile);
						 querySeqList2.add(seqPairRocket.reverse.recognizedSeqFile);
					  }
				   }
		       }
		   }
		
		   SeqMut SHM=new SeqMut();
		   SHM.setTmpDir(tmpDir);	
		   SHM.setAQS(minAlignLen,maxMismatch,maxGapNum,skipAQS);
	       SHM.setNQS(NQS_NN,NQS_minS,NQS_minNS,skipNQS);
	       SHM.setIndel(maxIndel,skipIndel);
	       SHM.setLMS(localBaseNum, localMaxMutNum, skipLMS);
	       SHM.setAvgRegion(avgRegion_start, avgRegion_end);
	       SHM.setSNP(snpRate, doSNPCheck);
	       SHM.setBaseErr(baseErrRate, doBaseErrCheck);
	       SHM.setTargetMut(interestMut);
		   SeqMutInfo mut;
		   List<SeqMutInfo> mutList;
		   String tmpFile=null;
		   if(seqType.equalsIgnoreCase(SeqOperation.SEQTYPE_SINGLEEND)){	
				if(querySeqList==null || querySeqList.size()==0){
					System.err.println("Error: No querySeq or querySeqList provided!");
					return;
				}
				mutList=new ArrayList<SeqMutInfo>();
				SHM.setIsReversed(false);
				for(int q=0;q<querySeqList.size();q++){
				   querySeqFile=querySeqList.get(q);			
				   if(SeqOperation.isFASTQSeq(querySeqFile) 
				    		  || SeqOperation.isFASTASeq(querySeqFile)){
					    
				    	outTag=querySeqFile.substring(
						   querySeqFile.lastIndexOf("/")+1,
						   (querySeqFile.lastIndexOf(".")>0) ? querySeqFile.lastIndexOf("."):querySeqFile.length()
					    );	
					    if(doSeqQCFilter && SeqOperation.isFASTQSeq(querySeqFile)){
						  System.out.println("Seq QC for ["+ querySeqFile+"]");	
						  seqQC=new SeqQCFilter();
						  seqQC.setOutDir(outDir);
						  if(seqQCOptList.size()>0){ 		
							  seqQCOpts=new String[seqQCOptList.size()+2];
							  for(int k=0;k<seqQCOptList.size();k++){
								 seqQCOpts[k]=seqQCOptList.get(k);
							  }
							  seqQCOpts[seqQCOptList.size()]="-out_format";
							  seqQCOpts[seqQCOptList.size()+1]="2";
						  }
						  seqQC.setOpts(seqQCOpts);
						  seqQC.prinseqQC(querySeqFile);
						  querySeqFile=seqQC.getResFasta();
						  seqQualList=seqQC.getSeqQual(seqQC.getResQual());	
						  tmpFile=querySeqFile;
						  if(seqQC.getResQual()!=null) FileOperate.delFile(seqQC.getResQual());						 
						  seqQC=null;
					    }else if(SeqOperation.isFASTQSeq(querySeqFile)){									 
			              querySeqFile=SeqOperation.convertFASTQ2FASTA(querySeqFile);
			              tmpFile=querySeqFile;
					    }					   
					
			            System.out.println("Call mutation for ["+ outTag+"]");			          
					    mut=SHM.getBLASTSeqMut(querySeqFile,refSeqFile,seqQualList);
					    seqQualList=null;
						if(tmpFile!=null) FileOperate.delFile(tmpFile);
						tmpFile=null;
					    if(mut!=null){ 						
						  SHM.saveMutInfo(mut, outDir, outTag);
						  mutList.add(mut);
						  mut=null;
					    }
					    if(!skipBaseMutCheck && doSNPCheck) SHM.saveSeqSNP(SHM.getSNPList(), outDir, outTag);
				   }else{
						System.out.println("Error: ["+querySeqFile+"] doesn't exist or isn't a fasta/fastq file. Skip mutation call for it!");
				   }				  
				}
				querySeqList=null;
				if(mutList.size()>0){
				  SHM.saveAvgBaseMutRate(mutList,outDir,"profile");
				  SHM.saveBaseMutRate(mutList,outDir,"profile");
				}
				mutList=null;
			
		   }else if(seqType.equalsIgnoreCase(SeqOperation.SEQTYPE_PAIREND)){				
				if(querySeqList==null || querySeqList.size()==0){
					System.err.println("Error: No querySeq or querySeqList provided!");
					return;
				}
				if(querySeqList2==null || querySeqList2.size()==0){
					System.err.println("Error: No querySeq2 or querySeqList2 provided!");
					return;
				}
			    String outTag2;
				SeqMutInfo mut2;
				List<SeqMutInfo> mutList2;
				mutList=new ArrayList<SeqMutInfo>();
				mutList2=new ArrayList<SeqMutInfo>();
					
				//forward
				SHM.setIsReversed(false);
				for(int q=0;q<querySeqList.size();q++){	
					  querySeqFile=querySeqList.get(q);
					  if(SeqOperation.isFASTQSeq(querySeqFile) 
					    		  || SeqOperation.isFASTASeq(querySeqFile)){
					     
						 outTag=querySeqFile.substring(
					       querySeqFile.lastIndexOf("/")+1,
					       (querySeqFile.lastIndexOf(".")>0) ? querySeqFile.lastIndexOf("."):querySeqFile.length()
					     );
						 
						 if(doSeqQCFilter && SeqOperation.isFASTQSeq(querySeqFile)){
							System.out.println("Seq QC for ["+ querySeqFile+"]");
							seqQC=new SeqQCFilter();
							seqQC.setOutDir(outDir);
							if(seqQCOptList.size()>0){ 		
							   seqQCOpts=new String[seqQCOptList.size()+2];
							   for(int k=0;k<seqQCOptList.size();k++){
								 seqQCOpts[k]=seqQCOptList.get(k);
							   }
							   seqQCOpts[seqQCOptList.size()]="-out_format";
							   seqQCOpts[seqQCOptList.size()+1]="2";
							}
							seqQC.setOpts(seqQCOpts);
							seqQC.prinseqQC(querySeqFile);
							querySeqFile=seqQC.getResFasta();
							seqQualList=seqQC.getSeqQual(seqQC.getResQual());
							tmpFile=querySeqFile;
						    if(seqQC.getResQual()!=null) FileOperate.delFile(seqQC.getResQual());
							seqQC=null;
						 }else if(SeqOperation.isFASTQSeq(querySeqFile)){									 
				            querySeqFile=SeqOperation.convertFASTQ2FASTA(querySeqFile);	
				            tmpFile=querySeqFile;
						 }
			
	                     System.out.println("Call mutation for ["+ outTag+"]");	                    
						 mut=SHM.getBLASTSeqMut(querySeqFile,refSeqFile,seqQualList);
						 seqQualList=null;
						 if(tmpFile!=null) FileOperate.delFile(tmpFile);
						 tmpFile=null;
						 if(mut!=null){ 						
						    SHM.saveMutInfo(mut, outDir, outTag);
						    mutList.add(mut);
						    mut=null;
					     }
						 if(!skipBaseMutCheck && doSNPCheck) SHM.saveSeqSNP(SHM.getSNPList(), outDir, outTag);
						 
					  }else{
						 System.out.println("Error: ["+querySeqFile+"] doesn't exist or isn't a fasta/fastq file. Skip mutation call for it!");
					  }					
				}
				querySeqList=null;
				if(mutList.size()>0){				
					SHM.saveBaseMutRate(mutList,outDir,"forward_profile");
				}
					  
				// reverse
				SHM.setIsReversed(true);
				for(int q=0;q<querySeqList2.size();q++){
					  querySeqFile2=querySeqList2.get(q);
					  if(SeqOperation.isFASTQSeq(querySeqFile2) 
					    		  || SeqOperation.isFASTASeq(querySeqFile2)){
						 
						 outTag2=querySeqFile2.substring(
						    querySeqFile2.lastIndexOf("/")+1,
						    (querySeqFile2.lastIndexOf(".")>0) ? querySeqFile2.lastIndexOf("."):querySeqFile2.length()
						 );	
						 
						 if(doSeqQCFilter && SeqOperation.isFASTQSeq(querySeqFile2)){											
							System.out.println("Seq QC for reverse: ["+ querySeqFile2+"]");
							seqQC=new SeqQCFilter();
							seqQC.setOutDir(outDir);
							if(seqQCOptList.size()>0){ 		
							   seqQCOpts=new String[seqQCOptList.size()+2];
							   for(int k=0;k<seqQCOptList.size();k++){
								 seqQCOpts[k]=seqQCOptList.get(k);
							   }
							   seqQCOpts[seqQCOptList.size()]="-out_format";
							   seqQCOpts[seqQCOptList.size()+1]="2";
							}
							seqQC.setOpts(seqQCOpts);
							seqQC.prinseqQC(querySeqFile2);
							querySeqFile2=seqQC.getResFasta();
							seqQualList2=seqQC.getSeqQual(seqQC.getResQual());
							tmpFile=querySeqFile2;
							if(seqQC.getResQual()!=null) FileOperate.delFile(seqQC.getResQual());
							seqQC=null;
						 }if(SeqOperation.isFASTQSeq(querySeqFile2)){									 
				            querySeqFile2=SeqOperation.convertFASTQ2FASTA(querySeqFile2);
				            tmpFile=querySeqFile2;
						 }
										
						 System.out.println("Call mutation for reverse: ["+ outTag2+"]");
						 mut2=SHM.getBLASTSeqMut(querySeqFile2,refSeqFile,seqQualList2);
						 seqQualList2=null;
						 if(tmpFile!=null) FileOperate.delFile(tmpFile);
						 tmpFile=null;
						 if(mut2!=null){ 						
						    SHM.saveMutInfo(mut2,outDir,outTag2);
						    mutList2.add(mut2);
						    mut2=null;
						 }
						 if(!skipBaseMutCheck && doSNPCheck) SHM.saveSeqSNP(SHM.getSNPList(),outDir,outTag);
						 
					   }else{
						 System.out.println("Error: ["+querySeqFile2+"] doesn't exist or isn't a fasta/fastq file. Skip mutation call for it!");
					   }					
				}//
				querySeqList2=null;
				if(mutList2.size()>0){
				  SHM.saveBaseMutRate(mutList2,outDir,"reverse_profile");
				}
				
				if(mutList.size()==mutList2.size() && mutList.size()>0){
			      SHM.saveAvgBaseMutRate(mutList,mutList2,outDir,"forward-reverse_profile");		
				}else{
				  if(mutList.size()>0){
					 SHM.saveAvgBaseMutRate(mutList,outDir,"forward_profile");		
				  }
				  if(mutList2.size()>0){
					 SHM.saveAvgBaseMutRate(mutList2,outDir,"reverse_profile");		
				  }
				}
				mutList2=null;
				
		   }//end if (single or pair end)		   
		   
		   //delTmpDir(tmpDir);	
		   if(tmpFiles!=null){
		    for(String tmp: tmpFiles){
			  FileOperate.delFile(tmp);
		    }
		   }
    }
}
