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
import org.geatools.data.structure.SeqInfo;
import org.geatools.operation.FileOperation;
import org.geatools.operation.SeqOperation;
import org.geatools.seqprocess.SeqFilter;
import org.geatools.seqprocess.SeqQCFilter;

public class CallUtility extends GEAT{
	
	static boolean doSplitSeq=false;	
	static int splitStep=1000000;
	static List<String> splitedSeqFiles = null;
	static List<String> splitedSeqFiles2 = null;
	
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
		 
		 String seqNameFile=null;
		 String seqNameCol=null;
		 int seqNameColIdx=0;
			
		 Map<String, List<String>> params=getCommandLineParams(args);
		 
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
		 }else if(params.get("fasta")!=null){
			isFastaOK=false;
			if(params.get("fasta").size()>0){
			  fasta=params.get("fasta").get(0).trim();
			  if(fasta!=null){
				 if(SeqOperation.isFASTASeq(fasta)){
				   isFastaOK=true;
				   fastaList=new ArrayList<String>();
				   fastaList.add(fasta);
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
		 
		 
		 if(params.get("fastq1")!=null){
			isFastqOK=false;
			if(params.get("fastq1").size()>0){
			   fastq=params.get("fastq1").get(0).trim();
			   if(fastq!=null){
			  	 if(SeqOperation.isFASTQSeq(fastq)){
			   	   isFastqOK=true;
			   	   fastqList=new ArrayList<String>();
			   	   fastqList.add(fastq);
			   	   seqType=SeqOperation.SEQTYPE_SINGLEEND;
			   	 }else{
			       System.err.println("The -fastq1 file doesn't exist or isn't a fastq file:(");
			       return;
			     }		    	
			   }
			}else{
			   System.err.println("You didn't provide fastq file :(");
			   return;
			}
		 }else if(params.get("fasta1")!=null){
			isFastaOK=false;
			if(params.get("fasta1").size()>0){
			  fasta=params.get("fasta1").get(0).trim();
			  if(fasta!=null){
				 if(SeqOperation.isFASTASeq(fasta)){
				   isFastaOK=true;
				   fastaList=new ArrayList<String>();
				   fastaList.add(fasta);
				   seqType=SeqOperation.SEQTYPE_SINGLEEND;
				 }else{
				   System.err.println("The -fasta1 file doesn't exist or isn't a fastq file:(");
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
			   	   fastaList2=new ArrayList<String>();
			   	   fastaList2.add(fasta2);
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
		 }else if(params.get("fastaList")!=null){
			isFastaOK=false;
			if(params.get("fastaList").size()>0){
				String fastqFiles=params.get("fastaList").get(0).trim();
				fastaList=FileOperation.getRowsOfFile(fastqFiles);
				if(SeqOperation.isFASTASeq(fastaList)){					
				   isFastaOK=true;
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
		 
		 
		 if(params.get("fastqList1")!=null){
			isFastqOK=false;
			if(params.get("fastqList1").size()>0){
			   String fastqFiles=params.get("fastqList1").get(0).trim();
		       fastqList=FileOperation.getRowsOfFile(fastqFiles);
		       if(SeqOperation.isFASTQSeq(fastqList)){					
				  isFastqOK=true;
			  	  seqType=SeqOperation.SEQTYPE_SINGLEEND;				
			   }else{	
				  System.err.println("Error: '-fastqList1' file doesn't exist or doesn't contain a fastq file.");
				  fastqList=new ArrayList<String>();					
			  	  return;
			   }
		  	   fastqFiles=null;
			}else{
			   System.err.println("Illegal '-fastqList1' parameter usage :(");			
			   return;
			}
		 }else if(params.get("fastaList1")!=null){
			isFastaOK=false;
			if(params.get("fastaList1").size()>0){
				String fastqFiles=params.get("fastaList1").get(0).trim();
				fastaList=FileOperation.getRowsOfFile(fastqFiles);
				if(SeqOperation.isFASTASeq(fastaList)){					
				   isFastaOK=true;
				   seqType=SeqOperation.SEQTYPE_SINGLEEND;				
				}else{	
				   System.err.println("Error: '-fastaList1' file doesn't exist or doesn't contain a fastq file.");
				   fastaList=new ArrayList<String>();					
				   return;
				}
				fastqFiles=null;
			}else{
				System.err.println("Illegal '-fastaList1' parameter usage :(");			
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
		 }else if(params.get("fastaList2")!=null){
			isFasta2OK=false;
			if(params.get("fastaList2").size()>0){
			    String fastqFiles=params.get("fastaList2").get(0).trim();					
			    fastaList2=FileOperation.getRowsOfFile(fastqFiles);
				if(SeqOperation.isFASTASeq(fastaList2)){
				   isFasta2OK=true;
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
		 
		 List<String> inSeqFiles = null;	
		 ////List<String> inSeqFiles2 = null;	
		 if(isFastaOK){
		    inSeqFiles=fastaList;
			////if(isFasta2OK) inSeqFiles2=fastaList2;
		 }else if(isFastqOK){
		    inSeqFiles=fastqList;
			////if(isFastq2OK) inSeqFiles2=fastqList2;
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
			   if(f.exists()){
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
		   
		   //####### count seq from given  fasta  
		 if(params.get("count")!=null){		
			  String seqFile = null;
			  if(isFastqOK || isFastaOK){
				System.out.println("Counting seq.........");
				if(isFastaOK){
				  for(int i=0;i<fastaList.size();i++){
				    seqFile=fastaList.get(i);
				    if(SeqOperation.isFASTASeq(seqFile)) 
				      System.out.println(SeqOperation.getSeqNum(seqFile) +" for ["+seqFile+"]");
					}
				}else if(isFastqOK){
				  for(int i=0;i<fastqList.size();i++){
					seqFile=fastqList.get(i);
					if(SeqOperation.isFASTQSeq(seqFile)) 
					  System.out.println(SeqOperation.getSeqNum(seqFile)+" for ["+seqFile+"]");
					}		
				 }
				 
			   }
			   
			   if(isFastq2OK || isFasta2OK){
				  System.out.println("Counting reverse seq.........");
				  if(isFasta2OK){
					for(int i=0;i<fastaList2.size();i++){
					    seqFile=fastaList2.get(i);
					    if(SeqOperation.isFASTASeq(seqFile)) 
					   	  System.out.println(SeqOperation.getSeqNum(seqFile));
						}
				  }else if(isFastq2OK){
					 for(int i=0;i<fastqList2.size();i++){
					   seqFile=fastqList2.get(i);
					   if(SeqOperation.isFASTQSeq(seqFile)) 
					   	 System.out.println(SeqOperation.getSeqNum(seqFile));
					 }		
				  }
			  }	
		 }
		   
		 //################Seq QC filter####################
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
				  System.out.println("Warning: empty '-seqQCFilter' parameter, just uses the default.");
				  doSeqQCFilter=true;
			   }				
		 } 
		   
		 //===================to do Seq QC filter=================
		 if(doSeqQCFilter && (isFastqOK || isFastq2OK)){
			  seqQC=new SeqQCFilter();
			  seqQC.setOutDir(outDir);
			  seqQC.setOpts(seqQCOpts);
			  String inSeqFile = fastq;	
			  if(seqType.equalsIgnoreCase(SeqOperation.SEQTYPE_SINGLEEND)){
				 if(SeqOperation.isFASTQSeq(inSeqFile)){ 
				   System.out.println("Seq QC for "+ inSeqFile);					    	
				   seqQC.prinseqQC(inSeqFile);
				   if(seqQC.getResFasta()!=null) fasta=seqQC.getResFasta();
				   else if(seqQC.getResFastq()!=null) fastq=seqQC.getResFastq();
				 }
			  }else if(seqType.equalsIgnoreCase(SeqOperation.SEQTYPE_PAIREND)){				
				 String inSeqFile2 = fastq2;
				 if(SeqOperation.isFASTQSeq(inSeqFile)){ 
				   System.out.println("Seq QC for "+ inSeqFile);						
				   seqQC.prinseqQC(inSeqFile);
				   if(seqQC.getResFasta()!=null) fasta=seqQC.getResFasta();
				   else if(seqQC.getResFastq()!=null) fastq=seqQC.getResFastq();
				 }
				 
				 if(SeqOperation.isFASTQSeq(inSeqFile2)){ 
				   System.out.println("Seq QC for reverse: "+ inSeqFile2);						
				   seqQC.prinseqQC(inSeqFile2);	
				   if(seqQC.getResFasta()!=null) fasta2=seqQC.getResFasta();
				   else if(seqQC.getResFastq()!=null) fastq2=seqQC.getResFastq();
				 }
			  }//SEQTYPE_PAIREND
		 }

		     
		 //####### split seq from given  fasta or fastq file########
		 //String splitedSeqOut=null;
		 if(params.get("split")!=null){
			doSplitSeq=true;				
			String str;			
			if(params.get("split").size()>0){
			   str=params.get("split").get(0).trim();
			   if(str.equalsIgnoreCase("no") || str.equalsIgnoreCase("n") 
					   || str.equalsIgnoreCase("false") || str.equalsIgnoreCase("f")) {
				  doSplitSeq=false;						 
			   }else{
				  if(NumberCheck.isPositiveInteger(str)){ 
					splitStep=Integer.parseInt(str); 
				  }else{
					System.out.println("Warning: '-split' value you set is illeagl. The default '-split "+splitStep +"' is forcefully used!");
					doSplitSeq=true;
				  }
			   }
			}else{
			   System.out.println("Warning: You didn't set '-split' value. The default '-split "+splitStep +"' is forcefully used!");
			   doSplitSeq=true;
		    }
		 }
		   
		 //=======================start to split seq===========================
		 if(doSplitSeq &&(isFastqOK || isFastaOK)){
			splitedSeqFiles=null;
			String inSeqFile="";
			//Check seq format, and then split seq into multiple subfiles................
			System.out.println("###### Spliting forward seq..........");
			if(isFastqOK){
			  inSeqFile=fastq;
			}else if(isFastaOK){	
			  inSeqFile=fasta;
		    }				
		
			if(outDir==null) outDir=inSeqFile.substring(0,inSeqFile.lastIndexOf("/"));
			String splitedSeqOut=outDir+"/split_forward";
			FileOperation.newFolder(splitedSeqOut);
	      	splitedSeqFiles=SeqOperation.splitSeqFile(inSeqFile,splitStep,splitedSeqOut);		 				
			if(splitedSeqFiles==null || splitedSeqFiles.size()==0){
			   System.err.println("Sorry, fail to split forward sequences!");
			   return;
			}
				
			if(seqType.equalsIgnoreCase(SeqOperation.SEQTYPE_PAIREND)){				
				splitedSeqFiles2=null;				
				String inSeqFile2="";
				//Check seq format, and then split seq into multiple subfiles................
				System.out.println("###### Spliting pair-end reverse seq.........");
				if(isFastq2OK && fastq2!=null){
				  inSeqFile2=fastq2;
				}else if(isFasta2OK && fasta2!=null){	
			      inSeqFile2=fasta2;
			    }		
				    
				if(outDir==null) outDir=inSeqFile2.substring(0,inSeqFile2.lastIndexOf("/"));
				splitedSeqOut=outDir+"/split_reverse";
				FileOperation.newFolder(splitedSeqOut);
			    splitedSeqFiles2=SeqOperation.splitSeqFile(inSeqFile2,splitStep,splitedSeqOut);				    	
			    if(splitedSeqFiles==null || splitedSeqFiles.size()==0){
			       System.err.println("Sorry, fail to split reverse sequences!");
				   return;
				}
			}				
		 }			   
		   
		 //####### extract seq by given seq name in file(s) ########
		 if(params.get("extractSeq")!=null){
			doSeqExtraction=false;
			String seqFile = null;	
			String format="fastq";
		    if(isFastaOK){
		       seqFile=fasta;
		       format="fasta";
		    }else if(isFastqOK){
		       seqFile=fastq;
		       format="fastq";
		    }
		    System.out.println("Total seq num: "+SeqOperation.getSeqNum(seqFile));
			if(outDir==null) outDir=seqFile.substring(0,seqFile.lastIndexOf("/"));
			
		    String outSeqFile;
		    if(outName==null) outName=seqFile.substring(				    
		    		   seqFile.lastIndexOf("/")+1,seqFile.lastIndexOf(".")
		    		 )+".extracted."+format;		    	
		    
		    outSeqFile=outDir+"/"+outName;
		    
			if(params.get("seqName") !=null && params.get("seqName").size()>1 
					 && (params.get("seqName").size() % 2 == 0)){
			   
			   System.out.println("================Extracting sub sequences from -seqName parameter===================");
			   List<SeqInfo> seqObjList = new ArrayList<SeqInfo>();
			   int k=0;
			   int seqNum0=0;		
			   for(int i=0;i<params.get("seqName").size();i=i+2){
				 k++;
				 seqNameFile=params.get("seqName").get(i).trim();
				 seqNameCol=params.get("seqName").get(i+1).trim();
			     if(seqNameFile!=null && NumberCheck.isPositiveInteger(seqNameCol)){
			    	doSeqExtraction=false;
			    	File f=new File(seqNameFile);
			    	if(f.exists()) doSeqExtraction=true;
			    	f=null;
			    	if(Integer.parseInt(seqNameCol)<=1) { 
			    	  seqNameColIdx=0;
			    	}else {
			    	  seqNameColIdx=Integer.parseInt(seqNameCol)-1;
			    	}
			    	System.out.println("Extracting for No."+k+" ......");
			    	if(doSeqExtraction){
			    	  SeqOperation.extractSubSeq(seqFile,seqNameFile,seqNameColIdx,outSeqFile);
			    	  seqObjList=SeqOperation.combineSeqList(seqObjList, SeqOperation.getSeqObj(outSeqFile));
			    	  System.out.println("Extracted seq: "+(seqObjList.size()-seqNum0));
			    	  System.out.println("+");
			    	  seqNum0=seqObjList.size();
			    	}else{
					  System.out.println("Warning: ["+seqNameFile+"] doesn't exist, we ignored it!!!");
					}		    	
			     }else{
					System.out.println("Warning: '-seqName' parameter is illeagl. We skiped it!!!");
				 }
			   }//for
			   if(seqObjList.size()>0 && k>1){
				 SeqOperation.saveSeqObj(seqObjList, outSeqFile);			  
			     System.out.println("Totally extracted seq: "+seqObjList.size());
			     System.out.println("Be Combined and saved in ["+outSeqFile+"]");
			     seqObjList=null;
			   }else if(seqObjList.size()==0){
				 System.err.println("No extracted sequences from ["+seqFile+"]");
				 return;
			   }
			}else{
			   System.out.println("Warning: '-seqName' parameter is illeagl. We skiped it!!!");
			}
		 }
		 
		 //####### extract seq name from given  fasta or fastq file########
		 if(params.get("extractSeqName")!=null){			 
			 System.out.println("================Extracting sequence name by parameter -extractSeqName===================");
			 boolean isOK=false;				
			 if(inSeqFiles!=null){				    		    	
			    String outFile;			    	  
				List<String> seqNameFiles = new ArrayList<String>();	
			    for(String seqFile:inSeqFiles){
			    	if(outDir==null) outDir=seqFile.substring(0,seqFile.lastIndexOf("/"));
					outFile=outDir+"/"+seqFile.substring(
					    		   seqFile.lastIndexOf("/")+1,seqFile.lastIndexOf(".")
					    		 )+".seqName";
			    		    		
			    	isOK=SeqOperation.extratSeqName(seqFile,outFile);
			    	if(isOK){
			    	   seqNameFiles.add(outFile);
			    	   System.out.println("Successfully extracted sequences name from ["+seqFile+"]");
			    	}else{
			    	   System.out.println("Error: failed to extract sequences name from ["+seqFile+"]");
			    	}
			    }
			    	  
			    if(outName!=null){
			    	String outCombined=outDir+"/"+outName;
					FileOperation.combineRowsOfFiles(seqNameFiles,outCombined);
					seqNameFiles=null;
					System.out.println("Combined sequences names and saved in ["+outCombined+"]");
					outCombined=null;
				}			  				    		    	
			    
			 }else {			    
		    	  System.out.println("Error: illegal parameter -extractSeqName");
		     }
			
		 }// -extractSeqName 

		   
		 //####### exclude seq by given seq name in file(s) ########
		 if(params.get("excludeSeq")!=null){
			doSeqExclusion=false;
			String inSeqFile = null;
			List<String> subSeqNameList=new ArrayList<String>();
			if(params.get("seqName") !=null && params.get("seqName").size()>1 
					 && (params.get("seqName").size() % 2 == 0)){
			   
			   System.out.println("================Excluding sub sequences according your -excludeSeq parameter===================");
			   List<SeqInfo> seqObjList = new ArrayList<SeqInfo>();		   
			   for(int i=0;i<params.get("seqName").size();i=i+2){
				 seqNameFile=params.get("seqName").get(i).trim();
				 seqNameCol=params.get("seqName").get(i+1).trim();
			     if(seqNameFile!=null && NumberCheck.isPositiveInteger(seqNameCol)){
			    	doSeqExclusion=false;
			    	File f=new File(seqNameFile);
			    	if(f.exists()) doSeqExclusion=true;
			    	f=null;
			    	if(Integer.parseInt(seqNameCol)<=1) {
			    	  seqNameColIdx=0;
			        }else {
			    	  seqNameColIdx=Integer.parseInt(seqNameCol)-1;
			        }
			     
			        if(isFastaOK){
			          inSeqFile=fasta;
			    	}else if(isFastqOK){
			    	  inSeqFile=fastq;
			    	}
			    	
			    	if(doSeqExclusion){		   
			    		subSeqNameList=SeqOperation.combineStrList(subSeqNameList,
			    				SeqOperation.getRowName(seqNameFile,seqNameColIdx));			    	 
			    	}		    	
			     }
			   }//for
			   if(subSeqNameList.size()>0){
				    seqObjList=SeqOperation.excludeSeq(inSeqFile,subSeqNameList);
				    subSeqNameList=null;
				 if(seqObjList.size()>0){
		    	    System.out.println("Successfully excluded "+seqObjList.size()+
		    	    		" sequences in ["+seqNameFile+"] from ["+inSeqFile+"]");
			     
			   	    if(outDir==null) outDir=inSeqFile.substring(0,inSeqFile.lastIndexOf("/"));
			        String outSeqFile=outDir+"/"+inSeqFile.substring(
			        		inSeqFile.lastIndexOf("/")+1,inSeqFile.lastIndexOf(".")
			    		 )+".excluded.fna";
			        SeqOperation.saveSeqObj(seqObjList, outSeqFile);
			        seqObjList=null;
			        isFastaOK=true;
			        fasta=outSeqFile;
				 }else{
				    System.err.println("No sequence exists for exclusion from "+inSeqFile);
				 }
			   }else{
				 System.err.println("No excluded sequences from "+inSeqFile);
				 return;
			   }
			}else{
			   System.out.println("Warning: '-seqName' parameter is illeagl. We skiped it!!!");
			}			
		 }
		 
		 if(params.get("excludeSeqAndDup")!=null){
			 boolean doExcludeSeq=false;
			 if(inSeqFiles==null) return;
			 for(String seqFile:inSeqFiles){
			   System.out.println("Total seq num: "+SeqOperation.getSeqNum(seqFile));
			   if(outDir==null) outDir=seqFile.substring(0,seqFile.lastIndexOf("/"));			     
			   if(params.get("excludeSeqAndDup").size()>1 && (params.get("excludeSeqAndDup").size() % 2 == 0)){
				   
				  System.out.println("================Removing sequences===================");				
				  int k=0;	
			      for(int i=0;i<params.get("excludeSeqAndDup").size();i=i+2){
					 k++;
					 seqNameFile=params.get("excludeSeqAndDup").get(i).trim();
					 seqNameCol=params.get("excludeSeqAndDup").get(i+1).trim();
				     if(seqNameFile!=null && NumberCheck.isPositiveInteger(seqNameCol)){
				    	doExcludeSeq=false;
				    	File f=new File(seqNameFile);
				    	if(f.exists()) doExcludeSeq=true;
				    	f=null;
				    	if(Integer.parseInt(seqNameCol)<=1) { 
				    	  seqNameColIdx=0;
				    	}else {
				    	  seqNameColIdx=Integer.parseInt(seqNameCol)-1;
				    	}
				    	System.out.println("Checking for No."+k+" ......");
				    	if(doExcludeSeq){				    	
				    	  System.out.println("Seq after removal: "
				    	         +SeqFilter.removeSeqAndItsDup(seqFile,seqNameFile,seqNameColIdx));				    			   
				    	}else{
						  System.out.println("Warning: ["+seqNameFile+"] doesn't exist, we ignored it!!!");
						}		    	
				     }else{
						System.out.println("Warning: '-excludeSeqAndDup' parameter is illeagl. We skiped it!!!");
					 }
				  }//for
		
			    }else{
				  System.out.println("Warning: '-excludeSeqAndDup' parameter is illeagl. We skiped it!!!");
			    }
			  }
				  
			  if(!doExcludeSeq) System.out.println("Warning: '-excludeSeqAndDup' parameter is illeagl. We skiped it!!!");
		 }

		 
		 //####### combined seq from given list of fasta or fastq files ########
		 if(params.get("combineSeq")!=null){
			 boolean doSeqCombine=false;		
			 String files=null;				 
			 List<String> seqFileList = null;
			 if(params.get("list") !=null && params.get("list").size()>0){
				files=params.get("list").get(0).trim();
				File f=new File(files);
	    	    if(f.exists()){
	    	      doSeqCombine=true;
	    	      seqFileList= FileOperation.getRowsOfFile(files);
	    	    }else{
	    	      System.out.println("Warning: listFile "+files+" doesn't exist.");
			      f=null;
				}				  
			 }else if(params.get("files") !=null && params.get("files").size()>0){				
				seqFileList=new ArrayList<String>();
				for(int i=0;i<params.get("files").size();i++){				
				  String file=params.get("files").get(i).trim();
	              File f=new File(file);
	    	      if(f.exists()){
	    	    	doSeqCombine=true;
	    	    	seqFileList.add(file);
	    	      }else {
	    	    	System.out.println("Warning:"+file+" doesn't exist.");
	    	      }
			      f=null;
			    }
			 }else{
				System.out.println("Warning: '-combineSeq' parameter is illeagl.");
				doSeqCombine=false;
			 }
			 
			 if(doSeqCombine){			   
				if(seqFileList!=null && seqFileList.size()>0){
				  String seqFile=seqFileList.get(0);
				  String format=FileOperation.getFileFormat(seqFile);
				  if(outDir==null) 
					 outDir=seqFile.substring(0,seqFile.lastIndexOf("/"))+"/combined";
				  FileOperation.newFolder(outDir);
				  if(outName==null) {										
					 outName=seqFile.substring(
							 seqFile.lastIndexOf("/")+1,seqFile.lastIndexOf(".")
					         )+"_combined."+format;					
				  }
				  
				  outName=outDir+"/"+outName;		

				  SeqOperation.combineSeqFile(seqFileList, outName);
				}
				seqFileList=null;
			 }
	     }		   
		  
	     //FileOperation.delFolder(tmpDir);
		 if(tmpFiles!=null){
		   for(String tmpFile: tmpFiles){
			  FileOperation.delFile(tmpFile);
		   }
		 }
		 delTmpDir(tmpDir);	
    }

}

