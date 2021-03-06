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
import org.geatools.data.structure.SeqInfo;
import org.geatools.operation.FileOperate;
import org.geatools.seqprocess.SeqOperation;
import org.geatools.seqprocess.SeqQCFilter;

public class CallUtility extends GEAT{
	
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
	
	public static void doWork(String[] args){		
		 
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
			isFastaOK=false;
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
		       fastqList=FileOperate.getRowListFromFile(fastqFiles);
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
				fastaList=FileOperate.getRowListFromFile(fastqFiles);
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
			   
		 if(params.get("fastqList2")!=null){
			isFastq2OK=false;
			if(params.get("fastqList2").size()>0){
			  String fastqFiles=params.get("fastqList2").get(0).trim();					
			  fastqList2=FileOperate.getRowListFromFile(fastqFiles);
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
			    fastaList2=FileOperate.getRowListFromFile(fastqFiles);
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
		 doSplitSeq=false;
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
		   
		 //=======================start to split seq===========================
		 if(doSplitSeq &&(isFastqOK || isFastaOK)){
				splitedSeqFiles=null;
				String inSeqFile="";
				//Check seq format, and then split seq into multiple subfiles................
				System.out.println(".........Spliting forward seq..........");
				if(isFastqOK){
					inSeqFile=fastq;
				}else if(isFastaOK){	
					inSeqFile=fasta;
		    	}				
				//System.out.println("Total Seq Num: "+SeqOperation.getSeqNum(inSeqFile));
	      	    splitedSeqFiles=SeqOperation.splitSeqFile(inSeqFile,splitStep,
	      	    		splitedSeqOut);
		 				
				if(splitedSeqFiles==null || splitedSeqFiles.size()==0){
					System.err.println("Sorry, not working for sequences split!");
				    return;
				}
				
				if(seqType.equalsIgnoreCase(SeqOperation.SEQTYPE_PAIREND)){				
					splitedSeqFiles2=null;				
					String inSeqFile2="";
					//Check seq format, and then split seq into multiple subfiles................
					System.out.println(".........Spliting pair-end reverse seq.........");
				    if(isFastq2OK && fastq2!=null){
					  inSeqFile2=fastq2;
					}else if(isFasta2OK && fasta2!=null){	
			    	  inSeqFile2=fasta2;
			    	}					
					//System.out.println("Total Seq Num(Pair-end reverse): "
			    	//                    +SeqOperation.getSeqNum(inSeqFile2)
			    	//                  );
			        splitedSeqFiles2=SeqOperation.splitSeqFile(inSeqFile2,splitStep,
			        		splitedSeqOut);			
			    	
			    	if(splitedSeqFiles==null || splitedSeqFiles.size()==0){
					   System.err.println("Sorry, not working for sequences split!");
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
		    String outSeqFile=outDir+"/"+seqFile.substring(
		    		   seqFile.lastIndexOf("/")+1,seqFile.lastIndexOf(".")
		    		 )+".extracted."+format;
		     
			if(params.get("seqName") !=null && params.get("seqName").size()>1 
					 && (params.get("seqName").size() % 2 == 0)){
			   
			   System.out.println("================Extracting sub sequences according your -seqName parameter===================");
			   List<SeqInfo> seqObjList = new ArrayList<SeqInfo>();
			   int k=0;
			   int seqNum0=0;		
			   for(int i=0;i<params.get("seqName").size();i=i+2){
				 k++;
				 seqNameFile=params.get("seqName").get(i).trim();
				 seqNameCol=params.get("seqName").get(i+1).trim();
			     if(seqNameFile!=null && isInteger(seqNameCol)){
			    	doSeqExtraction=false;
			    	File f=new File(seqNameFile);
			    	if(f.exists()) doSeqExtraction=true;
			    	f=null;
			    	if(Integer.parseInt(seqNameCol)<=1) 
			    	  seqNameColIdx=0;
			    	else
			    	  seqNameColIdx=Integer.parseInt(seqNameCol)-1;

			    	System.out.println("Extracting for No."+k+" ......");
			    	if(doSeqExtraction){
			    	  SeqOperation.extractSubSeq(seqFile,seqNameFile,0,outSeqFile);
			    	  seqObjList=SeqOperation.combineSeqList(
			    		 seqObjList, SeqOperation.getSeqObj(outSeqFile)
			    	  );
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
				 SeqOperation.saveSeqList(seqObjList, outSeqFile);			  
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
			     if(seqNameFile!=null && isInteger(seqNameCol)){
			    	doSeqExclusion=false;
			    	File f=new File(seqNameFile);
			    	if(f.exists()) doSeqExclusion=true;
			    	f=null;
			    	if(Integer.parseInt(seqNameCol)<=1) 
			    	  seqNameColIdx=0;
			    	else
			    	  seqNameColIdx=Integer.parseInt(seqNameCol)-1;
			    	
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
			        SeqOperation.saveSeqList(seqObjList, outSeqFile);
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
		 
		 //####### combined seq from given list of fasta or fastq files ########
		 if(params.get("combineSeq")!=null){
			 boolean doSeqCombine=false;		
			 String files=null;				 
			 List<String> seqFileList = null;
			 if(params.get("listFile") !=null && params.get("listFile").size()>0){
				files=params.get("seqFileList").get(0).trim();
				File f=new File(files);
	    	    if(f.exists()){
	    	      doSeqCombine=true;
	    	      seqFileList= FileOperate.getRowListFromFile(files);
	    	    }else{
	    	      System.out.println("Warning: listFile "+files+" doesn't exist.");
			      f=null;
				}				  
			 }else if(params.get("seqFiles") !=null && params.get("seqFiles").size()>0){				
				seqFileList=new ArrayList<String>();
				for(int i=0;i<params.get("seqFiles").size();i++){				
				  String file=params.get("seqFiles").get(i).trim();
	              File f=new File(file);
	    	      if(f.exists()){
	    	    	doSeqCombine=true;
	    	    	seqFileList.add(file);
	    	      }else System.out.println("Warning:"+file+" doesn't exist.");
			    	f=null;
			    }
			 }else{
				System.out.println("Warning: '-combineSeq' parameter is illeagl.");
				doSeqCombine=false;
			 }
			 
			 if(doSeqCombine){			   
				if(seqFileList!=null && seqFileList.size()>0){
				  String seqFile=seqFileList.get(0);
				  String format=FileOperate.getFileFormat(seqFile);
				  if(outDir==null) 
					 outDir=seqFile.substring(0,seqFile.lastIndexOf("/"))+"/combined";
				  FileOperate.newFolder(outDir);
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
		   
		 //####### extract seq name from given  fasta or fastq file########
		 if(params.get("extractSeqName")!=null){
			 String inSeqFiles=null;
			 List<String> seqFileList =new ArrayList<String>();
			 String outCombinedFile=null;			 		
			 
			 if(params.get("extractSeqName").size()>0){
				String subParams;
				String [] itemSplited;
				for(int i=0;i<params.get("extractSeqName").size();i++){
				  subParams=params.get("extractSeqName").get(i).trim();
				  itemSplited=subParams.split("=");
				  if(itemSplited.length>1){
					if(itemSplited[0].equalsIgnoreCase("seqFile")){
						inSeqFiles=itemSplited[1];
						File f=new File(inSeqFiles);
				    	if(!f.exists()){
				    	   System.err.println("Error: '-extractSeqName seqFile' doesn't exist.");
						   return;
				    	}
				    	f=null;
					}else if(itemSplited[0].equalsIgnoreCase("fileList")){
						inSeqFiles=itemSplited[1];
						File f=new File(inSeqFiles);
				    	if(!f.exists()){
				    	   System.err.println("Error: '-extractSeqName fileList' doesn't exist.");
						   return;
				    	}
				    	f=null;
						seqFileList=FileOperate.getRowListFromFile(inSeqFiles);
					}else if(itemSplited[0].equalsIgnoreCase("combinedOut")){
						outCombinedFile=itemSplited[1];
				    }
				  }else{
					System.err.println("Error: '-extractSeqName' parameter is illeagl.");
					return;
				  }	
				}
			 }else{
				System.err.println("Error: '-extractSeqName' parameter is illeagl.");
			    return;
			 }			  
			 
			 System.out.println("================Extracting sequence name by parameter -extractSeqName===================");
			 boolean isOK=false;				
			 if(inSeqFiles!=null || seqFileList.size()>0){	
			     if(seqFileList.size()>0) {			    	
			    	  String outFile;			    	  
					  List<String> seqNameFiles = new ArrayList<String>();	
			    	  for(String seqFile:seqFileList){
			    		if(outDir==null) outDir=seqFile.substring(0,seqFile.lastIndexOf("/"));
					    outFile=outDir+"/"+seqFile.substring(
					    		   seqFile.lastIndexOf("/")+1,seqFile.lastIndexOf(".")
					    		 )+".seqName";
			    		    		
			    		isOK=SeqOperation.extratSeqName(seqFile,outFile);
			    		if(isOK){
			    		   seqNameFiles.add(outFile);
			    		   System.out.println(
			    			 "Successfully extracted sequences name from ["+seqFile+"]");
			    		}else{
			    		   System.out.println(
			    			 "Error: failed to extract sequences name from ["+seqFile+"]");
			    		}
			    	  }
			    	  
			    	  if(outCombinedFile!=null){
						 FileOperate.combineRowListFromFiles(seqNameFiles,outCombinedFile);
						 seqNameFiles=null;
						 System.out.println("Successfully combined sequences name from ["
						    		    +inSeqFiles+"], and saved in ["+outCombinedFile+"]");
					  }	
			     }else if(inSeqFiles!=null){			    	
			    	  if(outDir==null) outDir=inSeqFiles.substring(0,inSeqFiles.lastIndexOf("/"));
			    	  String seqNameFile=outDir+"/"+inSeqFiles.substring(
							  inSeqFiles.lastIndexOf("/")+1,inSeqFiles.lastIndexOf(".")
					    		 )+".seqName";
			    	  isOK=SeqOperation.extratSeqName(inSeqFiles,seqNameFile);
			    	  if(isOK){			    		
			    		   System.out.println(
			    			 "Successfully extracted sequences name from ["+inSeqFiles+"]");
			    	  }else{
			    		   System.out.println(
			    			 "Error: failed to extract sequences name from ["+inSeqFiles+"]");
			    	  }
			     }  				    		    	
			    
			 }else {			    
		    	  System.out.println("Error: illegal parameter -extractSeqName");
		     }
			 seqFileList=null;
		 }// -extractSeqName 
		  
	     FileOperate.delFolder(tmpDir);
    }

}

