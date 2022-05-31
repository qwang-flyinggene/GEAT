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
import org.geatools.operation.FileOperation;
import org.geatools.operation.SeqOperation;
import org.geatools.seqprocess.SeqQCFilter;


public class CallFq2Fa extends GEAT{
	
	
	public static void doWork(String[] args) {	 
		 
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
			    if(!taskName.equalsIgnoreCase("Fq2Fa")){					  
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
		 }   

		 if(params.get("fastqList")!=null){
			 isFastqOK=false;
		     if(params.get("fastqList").size()>0){
				String fastqFiles=params.get("fastqList").get(0).trim();
				fastqList=FileOperation.getRowsOfFile(fastqFiles);
				if(SeqOperation.isFASTQSeq(fastqList)){					
					isFastqOK=true;
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
		 }
		   
		 if(params.get("fastqList2")!=null){
			  isFastq2OK=false;
		      if(params.get("fastqList2").size()>0){
		    	 String fastqFiles=params.get("fastqList2").get(0).trim();					
		    	 fastqList2=FileOperation.getRowsOfFile(fastqFiles);
				 if(SeqOperation.isFASTQSeq(fastqList2)){
					isFastq2OK=true;
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
			  outDir=homeDir+fileSeparator+"working";
			  dir=new File(outDir);
			  if(!dir.exists()) FileOperation.newFolder(outDir);	
		      dir=null;
			  doOutput=true;  			
		   }
		   
		   if(params.get("seqQCFilter")!=null){		
			  doSeqQCFilter=true;
			  String subParams;
			  String [] itemSplited;
			  if(params.get("seqQCFilter").size()>0){
				 subParams=params.get("seqQCFilter").get(0).trim();
				 if(subParams!=null && subParams.equalsIgnoreCase("skip"))			    		  
				    doSeqQCFilter=false;
			     else{	
				    List<String> seqQCOptList=new ArrayList<String>();
				    for(int i=0;i<params.get("seqQCFilter").size();i++){				
					  subParams=params.get("seqQCFilter").get(i).trim();
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
					  System.out.println("Warning: Illegal '-seqQCFilter' parameter");
					  System.out.println("Warning: The system uses default '-seqQCFilter' parameter.");
					  doSeqQCFilter=true;
				    }
			      }
			  }else{
				  System.out.println("Warning: empty '-seqQCFilter' parameter,the system uses default '-seqQCFilter' parameter.");
				  doSeqQCFilter=true;
			  }
				
		   }
		   
		  //===================Seq QC filter=================
		   if(doSeqQCFilter){
			  seqQC=new SeqQCFilter();
			  seqQC.setOutDir(outDir);
			  seqQC.setOpts(seqQCOpts);
			  String inSeqFile = null;
			  if(seqType.equalsIgnoreCase(SeqOperation.SEQTYPE_SINGLEEND)){
				 if(isFastqOK && fastqList.size()>0){
					for(String seqFile:fastqList){
					   inSeqFile=seqFile;				
					   if(SeqOperation.isFASTQSeq(inSeqFile)){
						  System.out.println("Seq QC for "+ inSeqFile);					    	
						  seqQC.prinseqQC(inSeqFile);
						  if(seqQC.getResFasta()==null && seqQC.getResFastq()!=null)
							 SeqOperation.convertFASTQ2FASTA(seqQC.getResFastq());
						  else{
							 SeqOperation.convertFASTQ2FASTA(inSeqFile); 
						  }
					   }	                  
				    }
				 }else{ 
					System.err.println("Error: missed/invalid '-fastq' parameter!!!");
				 }
			  }else if(seqType.equalsIgnoreCase(SeqOperation.SEQTYPE_PAIREND)){				
				 String inSeqFile2 = null;		
				 if((isFastqOK && fastqList!=null) || (isFastq2OK && fastqList2!=null)){
					
					if((isFastqOK && fastqList.size()>0)){
					   for(String seqFile:fastqList){
						inSeqFile=seqFile;												
						if(SeqOperation.isFASTQSeq(inSeqFile)){
							System.out.println("Seq QC for "+ inSeqFile);
							seqQC.prinseqQC(inSeqFile);	
							if(seqQC.getResFasta()==null && seqQC.getResFastq()!=null)
							  SeqOperation.convertFASTQ2FASTA(seqQC.getResFastq());
						    else{
							  SeqOperation.convertFASTQ2FASTA(inSeqFile); 
						    }
						}							
				      }// for 
					}
					
					if((isFastq2OK && fastqList2.size()>0)){
					  for(String seqFile:fastqList2){
						inSeqFile2=seqFile;												
						if(SeqOperation.isFASTQSeq(inSeqFile2)){										
							System.out.println("Seq QC for reverse: "+ inSeqFile2);
							seqQC.prinseqQC(inSeqFile2);
							if(seqQC.getResFasta()==null && seqQC.getResFastq()!=null)
							   SeqOperation.convertFASTQ2FASTA(seqQC.getResFastq());
							else{
							   SeqOperation.convertFASTQ2FASTA(inSeqFile2); 
							}
						}			
				      }// for 
					}
				 }else{ 
					System.err.println("Error: missed/invalid '-fastq or -fastq2' parameter!!!");
				 }
			  }//SEQTYPE_PAIREND
		   }else{
			  String inSeqFile = null;
			  if(seqType.equalsIgnoreCase(SeqOperation.SEQTYPE_SINGLEEND)){
				 if(isFastqOK && fastqList!=null && fastqList.size()>0){
					for(String seqFile:fastqList){
					   inSeqFile=seqFile;				
					   if(SeqOperation.isFASTQSeq(inSeqFile)){
						  SeqOperation.convertFASTQ2FASTA(inSeqFile);
						  System.out.println("OK for ["+inSeqFile+"]");
					   }else{
						  System.out.println("!!! Warning: ["+inSeqFile+"] is not a fastq file.");
				       }	                  
				    }
			     }else{ 
					System.err.println("Error: missed/invalid '-fastq' parameter!!!");
				 }
			  }else if(seqType.equalsIgnoreCase(SeqOperation.SEQTYPE_PAIREND)){				
				 String inSeqFile2 = null;		
				 if((isFastqOK && fastqList!=null) || (isFastq2OK && fastqList2!=null)){
						
					if((isFastqOK && fastqList.size()>0)){
					   for(String seqFile:fastqList){
						 inSeqFile=seqFile;												
						 if(SeqOperation.isFASTQSeq(inSeqFile)){
						    SeqOperation.convertFASTQ2FASTA(inSeqFile);
							System.out.println("OK for ["+inSeqFile+"]");
						 }else{
						    System.out.println("!!! Warning: ["+inSeqFile+"] is not a fastq file.");
					     }	 		
					   }// for 
					}
					if(isFastq2OK && fastqList2.size()>0){
					   for(String seqFile:fastqList2){
					     inSeqFile2=seqFile;												
						 if(SeqOperation.isFASTQSeq(inSeqFile2)){
						    SeqOperation.convertFASTQ2FASTA(inSeqFile2);
					        System.out.println("OK for ["+inSeqFile2+"]");
						 }else{
						    System.out.println("!!! Warning: ["+inSeqFile2+"] is not a fastq file.");
					     }	 		
					   }// for
					}
				 }else{ 
					System.err.println("Error: missed/invalid '-fastq or -fastq2' parameter!!!");
				 }
			  }//SEQTYPE_PAIREND 
		   }
		   		
		   if(tmpFiles!=null){
			 for(String tmpFile: tmpFiles){
				FileOperation.delFile(tmpFile);
			 }
		   }
		   delTmpDir(tmpDir);	
    }
}