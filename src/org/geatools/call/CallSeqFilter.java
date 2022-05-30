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
import java.util.ArrayList;

import java.util.List;
import java.util.Map;

import org.geatools.GEAT;
import org.geatools.NumberCheck;
import org.geatools.operation.FileOperation;
import org.geatools.operation.SeqOperation;
import org.geatools.seqprocess.SeqFilter;

public class CallSeqFilter extends GEAT{

	static boolean isDeDup=false;
	static int minSeqLen=-1;
	
	public static void doWork(String[] args){		
		 
		 if(homeDir==null) homeDir=GEAT.getHomeDir();			 
		 if(fileSeparator==null) fileSeparator = GEAT.getfileSeparator();		 
		 if(workingDir==null) workingDir=GEAT.getWorkingDir();	
		 File dir=null;
		 dir=new File(workingDir);
		 if(!dir.exists()) FileOperation.newFolder(workingDir);
		 
		 if(tmpDir==null) tmpDir=GEAT.getTmpDir();			
		 dir=new File(tmpDir);
		 if(!dir.exists()) FileOperation.newFolder(tmpDir);
		 dir=null;
			
		 Map<String, List<String>> params=getCommandLineParams(args);
		 if(params.get("task")!=null){	
			if(params.get("task").size()>0){
			  taskName=params.get("task").get(0).trim();
			  if(taskName!=null){
			    if(!taskName.equalsIgnoreCase("SeqFilter")){					  
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
		    	   seqType=SeqOperation.SEQTYPE_SINGLEEND;
		    	   fastqList=new ArrayList<String>();
		    	   fastqList.add(fastq);
		    	 }else{
			       System.err.println("The -fastq file doesn't exist or isn't a fastq file:(");
			       return;
			     }		    	
		       }
		    }else{
		       System.err.println("You didn't provide fastq file :(");
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
			    	  seqType=SeqOperation.SEQTYPE_PAIREND;
			    	  fastqList2=new ArrayList<String>();
			    	  fastqList2.add(fastq2);
			       }else{
				      System.err.println("The -fastq2 file doesn't exist or isn't a fastq file :(");
				      return;
				   }			    
			    }
			  }else{
			    System.err.println("You didn't provide fastq2 file :(");
			    return;
			  }
		   }   
		   
		   //####### set querySeq #######
		   if(params.get("fastqList")!=null){
			 isFastqOK=false;
		     if(params.get("fastqList").size()>0){
				String fastqFiles=params.get("fastqList").get(0).trim();
				fastqList=FileOperation.getRowsOfFile(fastqFiles);
				if(SeqOperation.isFASTQSeq(fastqList)){					
					seqType=SeqOperation.SEQTYPE_SINGLEEND;
					isFastqOK =true;					
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
		   }
		   
		   if(params.get("fastqList2")!=null){
			 isFastq2OK=false;
		     if(params.get("fastqList2").size()>0){
		    	 String fastqFiles=params.get("fastqList2").get(0).trim();					
		    	 fastqList2=FileOperation.getRowsOfFile(fastqFiles);
				 if(SeqOperation.isFASTQSeq(fastqList2)){					
					seqType=SeqOperation.SEQTYPE_PAIREND;
					isFastq2OK=true;				 
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
		   }
		   
		   
		   if(params.get("fasta")!=null){
			 isFastaOK=false;
			 if(params.get("fasta").size()>0){
			    fasta=params.get("fasta").get(0).trim();
				if(fasta!=null){
				  if(SeqOperation.isFASTASeq(fasta)){
				    isFastaOK =true;		
				    seqType=SeqOperation.SEQTYPE_SINGLEEND;
				    fastaList=new ArrayList<String>();
				    fastaList.add(fasta);
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
				   
			if(params.get("fasta2")!=null){
			   isFasta2OK=false;
		  	   if(params.get("fasta2").size()>0){
				  fasta2=params.get("fasta2").get(0).trim();
				  if(fasta2!=null){
				     if(SeqOperation.isFASTASeq(fasta2)){
				   	    isFasta2OK=true;
					    seqType=SeqOperation.SEQTYPE_PAIREND;
				   	    fastaList2=new ArrayList<String>();
					    fastaList2.add(fasta2);
					 }else{
					    System.err.println("The -fastq2 file doesn't exist or isn't a fastq file :(");
					    return;
					 }			    
				   }
			    }else{
				   System.err.println("You didn't provide fastq2 file :(");
				   return;
				}
		     }   
				   
				   //####### set querySeq #######
			 if(params.get("fastaList")!=null){	
				isFastaOK=false;
			    if(params.get("fastaList").size()>0){
				   String fastaFiles=params.get("fastaList").get(0).trim();
			       fastaList=FileOperation.getRowsOfFile(fastaFiles);
			   	   if(SeqOperation.isFASTASeq(fastaList)){					
					  seqType=SeqOperation.SEQTYPE_SINGLEEND;
					  isFastaOK=true;					
				   }else{	
					  System.err.println("Error: '-fastaList' file doesn't exist or doesn't contain a fastq file.");
				 	  fastaList=new ArrayList<String>();					
					  return;
				   }
				   fastaFiles=null;
				}else{
				   System.err.println("Illegal '-fastaList' parameter usage :(");			
				   return;
			    }
			  }
				   
			  if(params.get("fastaList2")!=null){		    	
				isFasta2OK=false;
				if(params.get("fastaList2").size()>0){
				   String fastaFiles=params.get("fastaList2").get(0).trim();					
				   fastaList2=FileOperation.getRowsOfFile(fastaFiles);
				   if(SeqOperation.isFASTASeq(fastaList2)){					
					  seqType=SeqOperation.SEQTYPE_PAIREND;
					  isFasta2OK=true;				 
				   }else{						
					  System.err.println("Error: '-fastaList2' file doesn't exist or doesn't contain a fastq file.");
				  	  fastaList2 = new ArrayList<String>();			
				      return;				    
				   }
				   fastaFiles=null;
				}else{
				   System.err.println("Illegal '-fastaList2' parameter usage :(");			
			  	   return;
				}
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
				  FileOperation.newFolder(outDir);
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
			 if(!dir.exists()) FileOperation.newFolder(outDir);	
		     dir=null;
			 doOutput=true;  		
		   }
		   
		   List<String> inSeqFiles = null;	
		   List<String> inSeqFiles2 = null;	
		   if(isFastaOK){
			 inSeqFiles=fastaList;
			 if(isFasta2OK) inSeqFiles2=fastaList2;
		   }else if(isFastqOK){
			 inSeqFiles=fastqList;
			 if(isFastq2OK) inSeqFiles2=fastqList2;
		   }
		   
		   //####### set isDeDup #######
		   if(params.get("deDup")!=null){			
			 isDeDup=true;		 
		   }
		   
		   //####### set minSeqLen #######
		   if(params.get("minLen")!=null){		  	
			 if(params.get("minLen").size()>0){
			    String str=params.get("minLen").get(0).trim();
				if(NumberCheck.isPositiveInteger(str)) {
				  minSeqLen=Integer.parseInt(str);
				}else{			   	
				  System.err.println("Illegal '-minLen' parameter usage, value of '-minLen' must be positive integer!");
				  return;
				}
			 }else{
				System.err.println("Illegal '-minLen' parameter usage!");
				return;
			 }
		   }
		   
		   ///===========================do it============================
		   if(inSeqFiles!=null && inSeqFiles2==null){
			 String seqFile;	
			 SeqFilter sFilter;
		     for(int i=0;i<inSeqFiles.size();i++){
		       seqFile=inSeqFiles.get(i).trim();
		       if(new File(seqFile).exists() 
		    	     && (SeqOperation.isFASTQSeq(seqFile) || SeqOperation.isFASTASeq(seqFile))) {
		    	  System.out.println("========Filtering sequences (file No."+(i+1)+")........");
		    	  sFilter=new SeqFilter();		    	 
		    	  sFilter.filter(seqFile,isDeDup, minSeqLen, outDir);
		    	  sFilter=null;
		       }else {
		    	  System.out.println("Warnning: No."+(i+1)+" ["+seqFile+ "] is illegal seq file, please check it!"); 
		       }
		     }
		   }if(inSeqFiles!=null && inSeqFiles2!=null && inSeqFiles.size()==inSeqFiles2.size()){
			 String seqFile;	
			 String seqFile2;
			 SeqFilter sFilter;
			 for(int i=0;i<inSeqFiles.size();i++){
			   seqFile=inSeqFiles.get(i).trim();
			   seqFile2=inSeqFiles2.get(i).trim();
			   if(new File(seqFile).exists() && new File(seqFile2).exists()) {
			   	  if((SeqOperation.isFASTQSeq(seqFile) && SeqOperation.isFASTQSeq(seqFile2))
			   		   || (SeqOperation.isFASTASeq(seqFile) && SeqOperation.isFASTASeq(seqFile2))
			   	    ){
			    	  System.out.println("========Filtering sequences (file No."+(i+1)+")........");
			    	  sFilter=new SeqFilter();			    	
			    	  sFilter.filter(seqFile,seqFile2,isDeDup, minSeqLen, outDir);
			    	  sFilter=null;
			   	  }
			   }else {
			   	 System.out.println("Warnning: No."+(i+1)+" illegal seq file, please check it!"); 
			   }
			 }  
		   }
		   
		   //delTmpDir(tmpDir);	
		   if(tmpFiles!=null){
			 for(String tmpFile: tmpFiles){
			   FileOperation.delFile(tmpFile);
			 }
		   }
    }

}
