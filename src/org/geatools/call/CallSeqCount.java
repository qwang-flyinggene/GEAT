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
import org.geatools.operation.FileOperate;
import org.geatools.seqprocess.SeqOperation;
import org.geatools.seqprocess.SeqQCFilter;
import org.geatools.seqprocess.QuerySeqCount;

public class CallSeqCount extends GEAT{
	
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
		 QuerySeqCount querySeqCount;
		 String querySeqInfo = null;
		 boolean isQuerySeqOK=false;
		 boolean skipNonExact=false;		 
		 splitStep=20000000;
		 doSplitSeq=false;
		 String inSeqFile=null;
		 String inSeqFile2=null;
		 List<String> inSeqFileList=null;	
		 List<String> inSeqFileList2=null;	
		 
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
			    if(!taskName.equalsIgnoreCase("SeqCount")){					  
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
			    	   fastqList=new ArrayList<String>();
			    	   fastqList.add(fastq);
			    	   inSeqFileList=fastqList;
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
				    	   fastaList=new ArrayList<String>();
				    	   fastaList.add(fasta);
				    	   inSeqFileList=fastaList;
				    	   seqType=SeqOperation.SEQTYPE_SINGLEEND;
				    	 }else{
					       System.err.println("The -fasta file doesn't exist or isn't a fastq file:(");
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
				    	   fastqList2=new ArrayList<String>();
				    	   fastqList2.add(fastq2);
				    	   inSeqFileList2=fastqList2;
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
				   isFastaOK=false;
				   if(params.get("fasta2").size()>0){
					   fasta2=params.get("fasta2").get(0).trim();
				       if(fasta2!=null){
				    	 if(SeqOperation.isFASTASeq(fasta2)){
				    	   isFasta2OK=true;
				    	   fastaList2=new ArrayList<String>();
				    	   fastaList2.add(fasta2);
				    	   inSeqFileList2=fastaList2;
				    	   seqType=SeqOperation.SEQTYPE_PAIREND;
				    	 }else{
					       System.err.println("The -fasta2 file doesn't exist or isn't a fastq file :(");
					       return;
					     }			    
				       }
					}else{
					   System.err.println("You didn't provide fasta2 file :(");
					   return;
					}
		  }   

		  if(params.get("fastqList")!=null){
				 isFastqOK=false;
			     if(params.get("fastqList").size()>0){
					String fastqFiles=params.get("fastqList").get(0).trim();
					fastqList=FileOperate.getRowListFromFile(fastqFiles);
					if(SeqOperation.isFASTQSeq(fastqList)){					
						isFastqOK=true;
						inSeqFileList=fastqList;
						seqType=SeqOperation.SEQTYPE_SINGLEEND;				
					}else{	
						System.err.println("Error: '-fastqList' file doesn't exist or doesn't contain a fastq file.");
						fastqList=new ArrayList<String>();					
						return;
					}
					fastqFiles=null;
				 }else{
					System.err.println("Illegal '-fastqList' parameter usage :(");			
					return;
				 }
		  }else if(params.get("fastaList")!=null){
				 isFastaOK=false;
			     if(params.get("fastaList").size()>0){
					String fastqFiles=params.get("fastaList").get(0).trim();
					fastaList=FileOperate.getRowListFromFile(fastqFiles);
					if(SeqOperation.isFASTASeq(fastaList)){					
						isFastaOK=true;
						inSeqFileList=fastaList;
						seqType=SeqOperation.SEQTYPE_SINGLEEND;				
					}else{	
						System.err.println("Error: '-fastaList' file doesn't exist or doesn't contain a fastq file.");
						fastaList=new ArrayList<String>();					
						return;
					}
					fastqFiles=null;
				 }else{
					System.err.println("Illegal '-fastaList' parameter usage :(");			
					return;
				 }
		   }
			   
		   if(params.get("fastqList2")!=null){
				  isFastq2OK=false;
			      if(params.get("fastqList2").size()>0){
			    	 String fastqFiles=params.get("fastqList2").get(0).trim();					
			    	 fastqList2=FileOperate.getRowListFromFile(fastqFiles);
					 if(SeqOperation.isFASTQSeq(fastqList2)){
						isFastq2OK=true;
						inSeqFileList2=fastqList2;
						seqType=SeqOperation.SEQTYPE_PAIREND;
					 }else{						
						System.err.println("Error: '-fastqList2' file doesn't exist or doesn't contain a fastq file.");
						fastqList2 = new ArrayList<String>();			
					    return;				    
					 }
					 fastqFiles=null;
				  }else{
					 System.err.println("Illegal '-fastqList2' parameter usage :(");			
					 return;
				  }
		  }else if(params.get("fastaList2")!=null){
				  isFasta2OK=false;
			      if(params.get("fastaList2").size()>0){
			    	 String fastqFiles=params.get("fastaList2").get(0).trim();					
			    	 fastaList2=FileOperate.getRowListFromFile(fastqFiles);
					 if(SeqOperation.isFASTASeq(fastaList2)){
						isFasta2OK=true;
						inSeqFileList2=fastaList2;
						seqType=SeqOperation.SEQTYPE_PAIREND;
					 }else{						
						System.err.println("Error: '-fastaList2' file doesn't exist or doesn't contain a fastq file.");
						fastaList2 = new ArrayList<String>();			
					    return;				    
					 }
					 fastqFiles=null;
				  }else{
					 System.err.println("Illegal '-fastaList2' parameter usage :(");			
					 return;
				  }
		  }		
		   
		   //####### set seq info for experiments in library#######
		  if(params.get("querySeq")!=null){
			  isQuerySeqOK=false;
			  if(params.get("querySeq").size()>0){
			     querySeqInfo=params.get("querySeq").get(0).trim();
			     if(querySeqInfo!=null){
			       File f = new File(querySeqInfo);
			       if(f.exists()){
			    	 isQuerySeqOK=true;
			    	 seqType=SeqOperation.SEQTYPE_SINGLEEND;
			       }else{
				     System.err.println(" '-querySeq' file you provided doesn't exist :(");
				     return;
				   }
			       f=null;
			     }
			  }else{
			     System.err.println("You didn't provide '-querySeq' file :(");
			     return;
			  }
		  }
		  		 
		  if(params.get("skip_nonexact")!=null){
			  skipNonExact=true;
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
		   
		  if(params.get("filter_seqQC")!=null){		
			   doSeqQCFilter=true;
			   String subParams;
			   String [] itemSplited;
			   if(params.get("filter_seqQC").size()>0){
				 subParams=params.get("filter_seqQC").get(0).trim();
				 if(subParams!=null && subParams.equalsIgnoreCase("skip"))			    		  
				    doSeqQCFilter=false;
			     else{	
				    List<String> seqQCOptList=new ArrayList<String>();
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
						}else if(itemSplited[0].trim().equalsIgnoreCase("out_format")){
							seqQCOptList.add("-out_format");
							seqQCOptList.add(itemSplited[1].trim());
						}else if(itemSplited[0].trim().equalsIgnoreCase("de_exact_dup")){
							seqQCOptList.add("-de_exact_dup");
							seqQCOptList.add(itemSplited[1].trim());
						}						
					  }
				    }
				    if(seqQCOptList.size()>0){ 
					  doSeqQCFilter=true;
					  seqQCOpts=new String[seqQCOptList.size()];
					  for(int k=0;k<seqQCOptList.size();k++){
						 seqQCOpts[k]=seqQCOptList.get(k);
					  }
				    }else{
					  System.out.println("Warning: Illegal '-filter_seqQC' parameter");
					  System.out.println("Warning: The system uses default '-filter_seqQC' parameter.");
					  doSeqQCFilter=true;
				    }
			      }
			   }else{
				  System.out.println("Warning: empty '-filter_seqQC' parameter,the system uses default '-filter_seqQC' parameter.");
				  doSeqQCFilter=true;
			   }				
		  }
		   
		  if(!isFastqOK){
			  System.out.println("Warning: no Fastq file provided, the system will skip seqQC step.");
			  doSeqQCFilter=false; 
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
		  
	
 		  if(isFastqOK && fastq!=null){
	    	   inSeqFile=fastq;
	      }else if(isFastaOK && fasta!=null){
	    	   inSeqFile=fasta;
	      } 	
 		 
 		  if(isFastq2OK && fastq2!=null){
	    	   inSeqFile2=fastq2;
	      }else if(isFasta2OK && fasta2!=null){
	    	   inSeqFile2=fasta2;
	      }
 		
		   //===================to do Seq QC filter=================
		  if(doSeqQCFilter && isFastqOK){
			  seqQC=new SeqQCFilter();
			  seqQC.setOutDir(outDir);
			  seqQC.setOpts(seqQCOpts);			
			  if(seqType.equalsIgnoreCase(SeqOperation.SEQTYPE_SINGLEEND)){
				 if(SeqOperation.isFASTQSeq(inSeqFile)){ 
				   System.out.println("Seq QC for ["+ inSeqFile+"]");	
				   System.out.println("Seq Counts: "+SeqOperation.getSeqNum(inSeqFile) +" for ["+inSeqFile+"]");
				   seqQC.prinseqQC(inSeqFile);
				   if(seqQC.getResFastq()!=null){ 
					   fastq=seqQC.getResFastq();
					   inSeqFile=fastq;
				   }else if(seqQC.getResFasta()!=null){
					   fasta=seqQC.getResFasta();
					   inSeqFile=fasta;
				   }
				   System.out.println("Seq Counts after QC: "+SeqOperation.getSeqNum(inSeqFile) +" for ["+inSeqFile+"]");					  
				 }
			  }else if(seqType.equalsIgnoreCase(SeqOperation.SEQTYPE_PAIREND)){					
				 if(SeqOperation.isFASTQSeq(inSeqFile)){
				   System.out.println("Seq QC for ["+ inSeqFile+"]");	
				   System.out.println("Seq Counts: "+SeqOperation.getSeqNum(inSeqFile) +" for ["+inSeqFile+"]");
				   seqQC.prinseqQC(inSeqFile);
				   if(seqQC.getResFastq()!=null){ 
					   fastq=seqQC.getResFastq();
					   inSeqFile=fastq;
				   }else if(seqQC.getResFasta()!=null){
					   fasta=seqQC.getResFasta();
					   inSeqFile=fasta;
				   }
				   System.out.println("Seq Counts after QC: "+SeqOperation.getSeqNum(inSeqFile) +" for ["+inSeqFile+"]");					  
				 }
				 
				 if(SeqOperation.isFASTQSeq(inSeqFile2)){ 
				   System.out.println("Seq QC for reverse: "+ inSeqFile2);	
				   System.out.println("Reverse Seq Counts: "+SeqOperation.getSeqNum(inSeqFile) +" for ["+inSeqFile+"]");
				   seqQC.prinseqQC(inSeqFile2);	
				   if(seqQC.getResFastq()!=null){ 
					   fastq2=seqQC.getResFastq();
					   inSeqFile2=fastq2;
				   }else if(seqQC.getResFasta()!=null){
					   fasta2=seqQC.getResFasta();
					   inSeqFile2=fasta2;
				   }
				   System.out.println("Reverse Seq Counts after QC: "+SeqOperation.getSeqNum(inSeqFile) +" for ["+inSeqFile+"]");					  
				 }
			  }//SEQTYPE_PAIREND
		  }
		   
		   //=======================start to split seq===========================
		  if(doSplitSeq &&(isFastqOK || isFastaOK)){
			  splitedSeqFiles=null;			
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
		   
		   //===================Recognizing reads for each querySeq=================
		  if((isFastqOK || isFastaOK) && isQuerySeqOK){
	    	 System.out.println("..........Recognizing reads for each querySeq..........");
	    	 isSeqRocketOK=false;
	    	 querySeqCount=new QuerySeqCount();    	   
	    	 querySeqCount.setHomeDir(homeDir);
	    	 querySeqCount.setDataDir(dataDir);
	    	 querySeqCount.setTmpDir(tmpDir);
	    	 if(skipNonExact) querySeqCount.setDoNonexact(false);
	    	 if(inSeqFileList!=null && inSeqFileList.size()>0){
		       if(seqType.equalsIgnoreCase(SeqOperation.SEQTYPE_SINGLEEND)){
		    	 querySeqCount.launchQueryCount(inSeqFileList,querySeqInfo,outDir,outTag,doSplitSeq);	    		
			     isSeqRocketOK=querySeqCount.isSeqRocketsOK();
		       }
		     }
		  }else{
			  //####### count seq from given  fasta /fastq 			  
			  if(isFastqOK || isFastaOK){
				System.out.println("Counting seq.........");
				for(int i=0;i<inSeqFileList.size();i++){
				  inSeqFile=inSeqFileList.get(i);
				  System.out.println(SeqOperation.getSeqNum(inSeqFile)+" for ["+inSeqFile+"]");
				}
			  }   
			  
			  if(isFastq2OK || isFasta2OK){
				System.out.println("Counting reverse seq.........");
			  	for(int i=0;i<inSeqFileList2.size();i++){
				  inSeqFile2=inSeqFileList2.get(i);
				  System.out.println(SeqOperation.getSeqNum(inSeqFile2)+" for ["+inSeqFile2+"]");
				}
			  }			  
		  }
		  
		  //delTmpDir(tmpDir);	
		  if(tmpFiles!=null){
			 for(String tmpFile: tmpFiles){
			   FileOperate.delFile(tmpFile);
			 }
		  }
    }

}
