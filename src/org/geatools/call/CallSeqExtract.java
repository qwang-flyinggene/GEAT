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
import org.geatools.seqprocess.SeqRocketConsole;

public class CallSeqExtract extends GEAT{
	
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
			    if(!taskName.equalsIgnoreCase("SeqExtract")){					  
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
		   
		   String inSeqFile = null;
	  	   String inSeqFile2 = null; // for reverse seq
	  	   
	  	   //if if both fastq and fasta are provided, we only analyze fastq.  
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
		   if(doSeqQCFilter && (isFastqOK || isFastq2OK)){
			  seqQC=new SeqQCFilter();
			  seqQC.setOutDir(outDir);
			  seqQC.setOpts(seqQCOpts);
			  if(seqType.equalsIgnoreCase(SeqOperation.SEQTYPE_SINGLEEND)){
				 if(SeqOperation.isFASTQSeq(inSeqFile)){ 
				   System.out.println("Seq QC for ["+ inSeqFile+"]");	
				   System.out.println("Seq Counts: "+SeqOperation.getSeqNum(inSeqFile) +" for ["+inSeqFile+"]");
				   seqQC.prinseqQC(inSeqFile);
				   if(seqQC.getResFasta()!=null){ 
					   fasta=seqQC.getResFasta();
					   inSeqFile=fasta;
				   } 
				   if(seqQC.getResFastq()!=null){
					   fastq=seqQC.getResFastq();
					   inSeqFile=fastq;
				   }
				   System.out.println("Seq Counts after QC: "+SeqOperation.getSeqNum(inSeqFile) +" for ["+inSeqFile+"]");					  
				 }
			  }else if(seqType.equalsIgnoreCase(SeqOperation.SEQTYPE_PAIREND)){	
				 if(SeqOperation.isFASTQSeq(inSeqFile)){
				   System.out.println("Seq QC for ["+ inSeqFile+"]");	
				   System.out.println("Seq Counts: "+SeqOperation.getSeqNum(inSeqFile) +" for ["+inSeqFile+"]");
				   seqQC.prinseqQC(inSeqFile);
				   if(seqQC.getResFasta()!=null){ 
					   fasta=seqQC.getResFasta();
					   inSeqFile=fasta;
				   }
				   if(seqQC.getResFastq()!=null){
					   fastq=seqQC.getResFastq();
					   inSeqFile=fastq;
				   }
				   System.out.println("Seq Counts after QC: "+SeqOperation.getSeqNum(inSeqFile) +" for ["+inSeqFile+"]");					  

				 }
				 
				 if(SeqOperation.isFASTQSeq(inSeqFile2)){ 
				   System.out.println("Seq QC for reverse: "+ inSeqFile2);	
				   System.out.println("Reverse Seq Counts: "+SeqOperation.getSeqNum(inSeqFile2) +" for ["+inSeqFile2+"]");
				   seqQC.prinseqQC(inSeqFile2);	
				   if(seqQC.getResFasta()!=null){ 
					   fasta2=seqQC.getResFasta();
					   inSeqFile2=fasta2;
				   }
				   if(seqQC.getResFastq()!=null){
					   fastq2=seqQC.getResFastq();
					   inSeqFile2=fastq2;
				   }
				   System.out.println("Reverse Seq Counts after QC: "+SeqOperation.getSeqNum(inSeqFile2) +" for ["+inSeqFile2+"]");					  

				 }
			  }//SEQTYPE_PAIREND
		   }
		   
		   //=======================start to split seq===========================
		   if(doSplitSeq &&(isFastqOK || isFastaOK)){
				//Check seq format, and then split seq into multiple subfiles................
				System.out.println(".........Working on forward seq..........");
				splitedSeqOut=outDir+"/split_forward";
				FileOperate.newFolder(splitedSeqOut);
			    splitedSeqFiles=null;			
	      	    splitedSeqFiles=SeqOperation.splitSeqFile(inSeqFile,splitStep,
	      	    		splitedSeqOut);
		 				
				if(splitedSeqFiles==null || splitedSeqFiles.size()==0){
					System.err.println("Sorry, not working for sequences split!");
				    return;
				}
				
				if(seqType.equalsIgnoreCase(SeqOperation.SEQTYPE_PAIREND)){				
					splitedSeqFiles2=null;
					//Check seq format, and then split seq into multiple subfiles................
					System.out.println(".........Working on pair-end reverse seq.........");
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

		    	}else if(seqType.equalsIgnoreCase(SeqOperation.SEQTYPE_PAIREND)){
		    		if(doSplitSeq)
		    		   seqRC.splitLaunchPairEnd(splitedSeqFiles,splitedSeqFiles2,
				        		expSeqInfo,expSeqInfo2,outDir);	
		    		else
		    		   seqRC.launchPairEnd(inSeqFile,inSeqFile2,expSeqInfo,expSeqInfo2,
		    					 outDir);
		    		   
			    	isSeqPairRocketOK=seqRC.isSeqRocketsOK();
			    }
		   }
	
		   if(tmpFiles!=null){
			 for(String tmpFile: tmpFiles){
			   FileOperate.delFile(tmpFile);
			 }
		   }
		   //delTmpDir(tmpDir);	
    }

}
