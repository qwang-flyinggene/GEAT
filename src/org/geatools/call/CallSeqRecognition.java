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
import org.geatools.seqprocess.SeqQCFilter;
import org.geatools.seqprocess.SeqRocketRecognition;
import org.geatools.seqprocess.SeqRocketConfig;
import org.geatools.seqprocess.SeqRocketConsole;

public class CallSeqRecognition extends GEAT{
	
	static boolean doSplitSeq=true;	
	static int splitStep=100000;
	static List<String> splitedSeqFiles = null;
	static List<String> splitedSeqFiles2 = null;
	
    static boolean finalSaveAsFASTA=false;
    static boolean finalSaveAsFASTQ=true;
    static boolean htmlSave=false;	    
	static String compoNameSaved=SeqRocketRecognition.BARCODE_NAME_DEFINITION;	
    static boolean compoSaveAsFASTA=false;
    static boolean compoSaveAsFASTQ=true; 
    static String seqCompoName5Trimmed=null;
    static String seqCompoName5Masked=null;
    static String seqCompoName3Trimmed=null;
    static String seqCompoName3Masked=null;
    static boolean mask5Save=false;	
    static boolean trim5Save=false;
    static boolean mask3Save=false;	
    static boolean trim3Save=false;
    
    static int barcodeMaxExactStart=1; 
    static float pctlMaxTStart3=0.75f;
    static int minAlignLen3=16;
    static boolean doFilterSeq=false; 
    static String FILTER_SEQTYPE_RAW="raw";
    static String FILTER_SEQTYPE_MASKED="masked";
    static String FILTER_SEQTYPE_TRIMMED="trimmed";   
    static String filterSeqType=FILTER_SEQTYPE_TRIMMED;
	static boolean doDeDup=false;
	static int minSeqLen=-1;
	
	static String libExpInfoPairMode=SeqRocketConsole.LIBSEQINFO_PEMODE_O2O;	
	//static List<SeqRocket>seqRockets;
	//static List<SeqRocketPair>seqPairRockets;
	
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
			    if(!taskName.equalsIgnoreCase("SeqRecognition")){					  
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
		 
		 
		 if(params.get("fastq1")!=null){
			 isFastqOK=false;
			 if(params.get("fastq1").size()>0){
		       fastq=params.get("fastq1").get(0).trim();
		       if(fastq!=null){
		    	 if(SeqOperation.isFASTQSeq(fastq)){
		    	   isFastqOK=true;
		    	   //doSeqQCFilter=true;
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
				    	doSeqQCFilter=false;
				    	seqType=SeqOperation.SEQTYPE_SINGLEEND;				    	
				     }else{
					    System.err.println("The -fasta1 file doesn't exist or isn't a fasta file:(");
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
		   if(params.get("expInfo")!=null){
			  isLibExpSeqInfoOK=false;
			  if(params.get("expInfo").size()>0){
			     expSeqInfo=params.get("expInfo").get(0).trim();
			     if(expSeqInfo!=null){
			       File f = new File(expSeqInfo);
			       if(f.exists()){
			    	 isLibExpSeqInfoOK=true;
			    	 seqType=SeqOperation.SEQTYPE_SINGLEEND;
			       }else{
				     System.err.println(" '-expInfo' file you provided doesn't exist :(");
				     return;
				   }
			       f=null;
			     }
			  }else{
			     System.err.println("You didn't provide '-expInfo' file :(");
			     return;
			  }
		   }
		   
		   if(params.get("expInfo1")!=null){
			  isLibExpSeqInfoOK=false;
			  if(params.get("expInfo1").size()>0){
			     expSeqInfo=params.get("expInfo1").get(0).trim();
			     if(expSeqInfo!=null){
			       File f = new File(expSeqInfo);
			       if(f.exists()){
			    	 isLibExpSeqInfoOK=true;
			    	 seqType=SeqOperation.SEQTYPE_SINGLEEND;
			       }else{
				     System.err.println(" '-expInfo1' file you provided doesn't exist :(");
				     return;
				   }
			       f=null;
			     }
			  }else{
			     System.err.println("You didn't provide '-expInfo1' file :(");
			     return;
			  }
		   }

		   
		   if(params.get("expInfo2")!=null){
			  isLibExpSeqInfoOK=false;
			  if(params.get("expInfo2").size()>0){
			     expSeqInfo2=params.get("expInfo2").get(0).trim();
			     if(expSeqInfo2!=null){
			       File f = new File(expSeqInfo2);
			       if(f.exists()){
			    	 isLibExpSeqInfoOK=true;
			    	 seqType=SeqOperation.SEQTYPE_PAIREND;
			       }else{
				     System.err.println(" '-expInfo2' file you provided doesn't exist :(");
				     return;
				   }
			       f=null;
			     }
			  }else{
			     System.err.println("You didn't provide '-expInfo2' file :(");
			     return;
			  }
		   }		     
		   
		   // '-expInfoPEMode' value: 'pw' denotes 'pairwise' or 'o2o' denotes 'one2one'
		   if(params.get("expInfo_PEMode")!=null){			
			  if(params.get("expInfo_PEMode").size()>0){
				 libExpInfoPairMode=params.get("expInfo_PEMode").get(0).trim();
				 if(libExpInfoPairMode.equalsIgnoreCase(SeqRocketConsole.LIBSEQINFO_PEMODE_PW)
					|| libExpInfoPairMode.equalsIgnoreCase(SeqRocketConsole.LIBSEQINFO_PEMODE_O2O)) {
				     
					 System.out.println("'-expInfo_PEMode' "+libExpInfoPairMode);
				 }else{
					 System.err.println("Illegal '-expInfo_PEMode' parameter value, you must set it as 'pw' for 'pairwise' or o2o for 'one to one paired'.");
					 return;
				 }
			  }else{
				 System.err.println("Illegal '-expInfo_PEMode' parameter value, you must set it as 'pw' for 'pairwise' or o2o for 'one to one paired'.");
				 return;
			  }
		   }

		   
		   ////####### set output #######
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

		   
		   if(params.get("trim5")!=null){
			  trim5Save=true;
			  seqCompoName5Trimmed=null;
			  if(params.get("trim5").size()>0){
				 seqCompoName5Trimmed=params.get("trim5").get(0).trim();
				 if(seqCompoName5Trimmed.equalsIgnoreCase("no") 
						 || seqCompoName5Trimmed.equalsIgnoreCase("none") ) {
					trim5Save=false; 
					seqCompoName5Trimmed=null;
				 }
			  }
		   }
		   
		   if(params.get("mask5")!=null){
			  mask5Save=true;
			  seqCompoName5Masked=null;
			  if(params.get("mask5").size()>0){
				 seqCompoName5Masked=params.get("mask5").get(0).trim();
				 if(seqCompoName5Masked.equalsIgnoreCase("no") 
						 || seqCompoName5Masked.equalsIgnoreCase("none") ) {
					mask5Save=false; 
					seqCompoName5Masked=null;
				 }
			  }
		   }
		   
		   if(params.get("trim3")!=null){
			  trim3Save=true;
			  seqCompoName3Trimmed=null;
			  if(params.get("trim3").size()>0){
				 seqCompoName3Trimmed=params.get("trim3").get(0).trim();
				 if(seqCompoName3Trimmed.equalsIgnoreCase("no") 
						 || seqCompoName3Trimmed.equalsIgnoreCase("none") ) {
					trim3Save=false; 
					seqCompoName3Trimmed=null;
				 }
			  }
		   }
		   
		   if(params.get("mask3")!=null){
			  mask3Save=true;
			  seqCompoName3Masked=null;
			  if(params.get("mask3").size()>0){
				 seqCompoName3Masked=params.get("mask3").get(0).trim();
				 if(seqCompoName3Masked.equalsIgnoreCase("no") 
						 || seqCompoName3Masked.equalsIgnoreCase("none") ) {
					 mask3Save=false; 
					 seqCompoName3Masked=null;
				 }
			  }
		   }
		   
		   if(params.get("5_max_sstart_exact")!=null){		  	
			  if(params.get("5_max_sstart_exact").size()>0){
				  String str=params.get("5_max_sstart_exact").get(0).trim();
				  if(NumberCheck.isPositiveInteger(str)) {
					 barcodeMaxExactStart=Integer.parseInt(str);
			      }else {
		    		 System.err.println("'-5_max_sstart_exact' value must be a positive integer!!!");
		    		 return;
		          }
			  }else{
				  System.err.println("Illegal '-5_max_sstart_exact' parameter usage :(");
				  return;
			  }
		   }
		   
		   if(params.get("3_max_tstart_p")!=null){		  	
			  if(params.get("3_max_tstart_p").size()>0){
				  String str=params.get("3_max_tstart_p").get(0).trim();
				  if(NumberCheck.isPercentile(str)){
					 pctlMaxTStart3=Float.parseFloat(str);
			      }else {
		    		 System.err.println("'-3_max_tstart_p' value must be between 0 and 1!!!");
		    		 return;
		          }
			  }else{
				  System.err.println("Illegal '-3_max_tstart_p' parameter usage :(");
				  return;
			  }
		   }
		   
		   if(params.get("3_min_alignLen")!=null){		  	
			  if(params.get("3_min_alignLen").size()>0){
				  String str=params.get("3_min_alignLen").get(0).trim();
				  if(NumberCheck.isPositiveNumeric(str)){
					 minAlignLen3=Integer.parseInt(str);
			      }else {
		    		 System.err.println("'-3_min_alignLen' value must be a positive integer!!!");
		    		 return;
		          }
			  }else{
				  System.err.println("Illegal '-3_min_alignLen' parameter usage :(");
				  return;
			  }
		   }

		   
		   if(params.get("compoSave")!=null){		  	
			   if(params.get("compoSave").size()==2){				    	  
				   compoNameSaved=params.get("compoSave").get(0).trim();		
			       String str=params.get("compoSave").get(1).trim();
			       if(str.equalsIgnoreCase("fastq") || str.equalsIgnoreCase("fq")){
			    	  compoSaveAsFASTQ=true;
			    	  compoSaveAsFASTA=false;
			       }else if(str.equalsIgnoreCase("fasta") || str.equalsIgnoreCase("fa")) {			    		
			    	  compoSaveAsFASTQ=false;
			    	  compoSaveAsFASTA=true;
			       }else if(str.equalsIgnoreCase("none") || str.equalsIgnoreCase("no")){
			    	  compoSaveAsFASTA=false;
			    	  compoSaveAsFASTQ=false;
			       }else {
			    	  System.err.println("Illegal '-compoSave' parameter usage, we just use the default value");
			       }					 	   					 
				}else {
				   System.err.println("Illegal '-compoSave' parameter usage, we just use the default value");
				   return;
				}
		   }
		   
		   if(params.get("finalSave")!=null){		  	
			    if(params.get("finalSave").size()>0){			    
			      String str="";
			      for(int i=0;i<params.get("finalSave").size();i++) {
			    	  
			    	 str=params.get("finalSave").get(i).trim();			    	 

			    	 if(str.equalsIgnoreCase("fastq") || str.equalsIgnoreCase("fq")){
			    		 finalSaveAsFASTQ=true;
			    		 finalSaveAsFASTA=false;
			    	 }else if(str.equalsIgnoreCase("fasta") || str.equalsIgnoreCase("fa")) {
			    		 finalSaveAsFASTQ=false;
			    		 finalSaveAsFASTA=true;			    		 
			    	 }else if(str.equalsIgnoreCase("noSave")){
			    		 finalSaveAsFASTQ=false;
			    		 finalSaveAsFASTA=false;
			    	 }		    	 
			    	 
			    	 if(str.equalsIgnoreCase("html")){
			    		 htmlSave=true;
			    	 }else if(str.equalsIgnoreCase("noHTML")){
			    		 htmlSave=false;	
			    	 }

			      }			 	   					 
				}else{
				  System.err.println("Illegal '-finalSave' parameter usage :(");
				  return;
				}
		   }

		   
		   if(params.get("seqQCFilter")!=null){		
			   doSeqQCFilter=true;
			   String subParams;
			   String [] itemSplited;
			   if(params.get("seqQCFilter").size()>0){
				 subParams=params.get("seqQCFilter").get(0).trim();
				 if(subParams!=null && subParams.equalsIgnoreCase("skip")){		    		  
				    doSeqQCFilter=false;
			     }else{	
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
		   
		   if(!isFastqOK){
			  System.out.println("Warning: no Fastq file provided, the system will skip seqQC step.");
			  doSeqQCFilter=false; 
		   }
		   
		   //####### split seq from given  fasta or fastq file########
		   //String splitedSeqOut=null;
		   if(params.get("split")!=null){
			  doSplitSeq=true;				
			  String str;			
			  if(params.get("split").size()>0){
				 str=params.get("split").get(0).trim();
				 if(str.equalsIgnoreCase("no") || str.equalsIgnoreCase("n")) doSplitSeq=false;						 
				 else{
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
		   
		   //####### set isFilterSeq #######
		   if(params.get("filter")!=null){			
			 doFilterSeq=true;
			 if(params.get("filter").size()>0){
			   String str=params.get("filter").get(0).trim();
			   if(str.equalsIgnoreCase("no") || str.equalsIgnoreCase("n")) { 
				  doFilterSeq=false;	
			   }else if(str.equalsIgnoreCase(FILTER_SEQTYPE_RAW)) { //"raw"
				  filterSeqType=FILTER_SEQTYPE_RAW;
			   }else if(str.equalsIgnoreCase(FILTER_SEQTYPE_TRIMMED)) { //"trimmed"
				  filterSeqType=FILTER_SEQTYPE_TRIMMED;
			   }else if(str.equalsIgnoreCase(FILTER_SEQTYPE_MASKED)) { //"masked"
				  filterSeqType=FILTER_SEQTYPE_MASKED;
			   }
			 }
		   }

		   //####### set isDeDup #######
		   if(params.get("deDup")!=null){			
			 doDeDup=true;
			 if(params.get("deDup").size()>0){
			   String str=params.get("deDup").get(0).trim();
			   if(str.equalsIgnoreCase("no") || str.equalsIgnoreCase("n")) doDeDup=false;	
			 }
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
			  System.out.println("###### Spliting forward seq..........");
			  String splitedSeqOut=outDir+"/split_forward";
			  FileOperation.newFolder(splitedSeqOut);
			  splitedSeqFiles=null;			
	      	  splitedSeqFiles=SeqOperation.splitSeqFile(inSeqFile,splitStep,splitedSeqOut);		 				
			  if(splitedSeqFiles==null || splitedSeqFiles.size()==0){
				 System.err.println("Sorry, fail to split forward sequences!");
				 return;
			  }
				
			  if(seqType.equalsIgnoreCase(SeqOperation.SEQTYPE_PAIREND)){					
				  //Check seq format, and then split seq into multiple subfiles................
				 System.out.println("###### Spliting pair-end reverse seq.........");
				 splitedSeqOut=outDir+"/split_reverse";
				 FileOperation.newFolder(splitedSeqOut);
				 splitedSeqFiles2=null;
			     splitedSeqFiles2=SeqOperation.splitSeqFile(inSeqFile2,splitStep,splitedSeqOut);			    	
			     if(splitedSeqFiles2==null || splitedSeqFiles2.size()==0){
				    System.err.println("Sorry, fail to split reverse sequences!");
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
	    	   seqRC.setTmpFiles(tmpFiles);	
	    	   SeqRocketConfig rocketConfig=new SeqRocketConfig();
	    	   rocketConfig.setTmpDir(tmpDir);
	    	   rocketConfig.setTmpFiles(tmpFiles);	
	    	   rocketConfig.setBarcodeMaxExactStart(barcodeMaxExactStart);
	    	   rocketConfig.setRC3MaxTStartPCTL(pctlMaxTStart3);
	    	   rocketConfig.setRC3MinAlignLen(minAlignLen3);
	    	   rocketConfig.setTrimLeft(seqCompoName5Trimmed,trim5Save);
	    	   rocketConfig.setMaskLeft(seqCompoName5Masked,mask5Save);
	    	   rocketConfig.setTrimRight(seqCompoName3Trimmed,trim3Save);
	    	   rocketConfig.setMaskRight(seqCompoName3Masked,mask3Save);	    	  
	    	   rocketConfig.setSaveRecognizedSeqFormat(finalSaveAsFASTA,finalSaveAsFASTQ);
	    	   rocketConfig.setSaveHTML(htmlSave);
	    	   rocketConfig.setSaveSeqComponent(compoNameSaved,compoSaveAsFASTA,compoSaveAsFASTQ);
	    	   
	    	   if(seqType.equalsIgnoreCase(SeqOperation.SEQTYPE_SINGLEEND)){
		    	  if(doSplitSeq) {
		    		seqRC.splitLaunchSingleEnd(splitedSeqFiles,expSeqInfo,rocketConfig,outDir);
		    	  }else {
		    		seqRC.launchSingleEnd(inSeqFile,expSeqInfo,rocketConfig,outDir);
		    	  }	 
				  isSeqRocketOK=seqRC.isSeqRocketsOK();
				  
				  if(doFilterSeq){
					List<String> seqFileList = null;
					if(filterSeqType.equalsIgnoreCase(FILTER_SEQTYPE_RAW)) { // "raw"
					  seqFileList=seqRC.getRecognizedFileList();
					}else if(filterSeqType.equalsIgnoreCase(FILTER_SEQTYPE_TRIMMED)) { //"trimmed"
					  seqFileList=seqRC.getRecognizedTrimFileList();
					}if(filterSeqType.equalsIgnoreCase(FILTER_SEQTYPE_MASKED)) { //"masked"
					  seqFileList=seqRC.getRecognizedMaskFileList();
					}
					if(seqFileList!=null) {
				      for(int i=0;i<seqFileList.size();i++){
					    String seqFile=seqFileList.get(i).trim();
					    if(new File(seqFile).exists() 
					         && (SeqOperation.isFASTQSeq(seqFile) || SeqOperation.isFASTASeq(seqFile))) {
					   	  System.out.println("========Filtering sequences (file No."+(i+1)+")........");
				    	  SeqFilter sf=new SeqFilter();
				    	  if(doSplitSeq) sf.setSplitStep(splitStep);
					   	  sf.filter(seqFile,doDeDup, minSeqLen, outDir);
					   	  sf=null;
					    }else {
					   	  System.out.println("Warnning: No."+(i+1)+" ["+seqFile+ "] is illegal seq file, please check it!"); 
					    }
					  }
					}
					seqFileList=null;
				  }
		       }else if(seqType.equalsIgnoreCase(SeqOperation.SEQTYPE_PAIREND)){		    	  
			      SeqRocketConfig rocketConfig2=new SeqRocketConfig();
		    	  rocketConfig2.setTmpDir(tmpDir);
		    	  rocketConfig2.setTmpFiles(tmpFiles);	
		    	  rocketConfig2.setBarcodeMaxExactStart(barcodeMaxExactStart);
		    	  rocketConfig2.setRC3MaxTStartPCTL(pctlMaxTStart3);
		    	  rocketConfig2.setRC3MinAlignLen(minAlignLen3);
		    	  rocketConfig2.setTrimLeft(seqCompoName5Trimmed,trim5Save);
		    	  rocketConfig2.setMaskLeft(seqCompoName5Masked,mask5Save);
		    	  rocketConfig2.setTrimRight(seqCompoName3Trimmed,trim3Save);
		    	  rocketConfig2.setMaskRight(seqCompoName3Masked,mask3Save);	    	  
		    	  rocketConfig2.setSaveRecognizedSeqFormat(finalSaveAsFASTA,finalSaveAsFASTQ);
		    	  rocketConfig2.setSaveHTML(htmlSave);
		    	  rocketConfig2.setSaveSeqComponent(compoNameSaved,compoSaveAsFASTA,compoSaveAsFASTQ);
		    	
			      if(doSplitSeq){
			    	 seqRC.splitLaunchPairEnd(splitedSeqFiles,splitedSeqFiles2,
			    			  expSeqInfo,expSeqInfo2,libExpInfoPairMode,rocketConfig,rocketConfig2,outDir);	
		    	  }else{		    		
		    		 seqRC.launchPairEnd(inSeqFile,inSeqFile2,libExpInfoPairMode,
		    			     expSeqInfo,expSeqInfo2,rocketConfig,rocketConfig2,outDir);		    		 
		    	  }   
			      isSeqPairRocketOK=seqRC.isSeqRocketsOK();
			      
				  if(doFilterSeq) {
					List<ArrayList<String>> seqPairFileList = null;
					List<String> seqFileList1 = null;
					List<String> seqFileList2 = null;
					if(filterSeqType.equalsIgnoreCase(FILTER_SEQTYPE_RAW)) { //"raw"
					   seqPairFileList=seqRC.getRecognizedPairFileList();
					   if(seqPairFileList!=null) {
						 seqFileList1=seqPairFileList.get(0);
						 seqFileList2=seqPairFileList.get(1);
					   }
					}else if(filterSeqType.equalsIgnoreCase(FILTER_SEQTYPE_TRIMMED)) { //"trimmed"
					   seqPairFileList=seqRC.getRecognizedPairTrimFileList();
					   if(seqPairFileList!=null) {
						 seqFileList1=seqPairFileList.get(0);
						 seqFileList2=seqPairFileList.get(1);
					   }
					}if(filterSeqType.equalsIgnoreCase(FILTER_SEQTYPE_MASKED)) { //"masked"
					   seqPairFileList=seqRC.getRecognizedPairMaskFileList();
					   if(seqPairFileList!=null) {
						 seqFileList1=seqPairFileList.get(0);
						 seqFileList2=seqPairFileList.get(1);
					   }
					}
					
					if(seqPairFileList!=null && seqFileList1.size()==seqFileList2.size()){
					  for(int i=0;i<seqFileList1.size();i++){
					    String seqFile1=seqFileList1.get(i).trim();
					    String seqFile2=seqFileList2.get(i).trim();
					    if(new File(seqFile1).exists() && new File(seqFile2).exists()) {
						  if((SeqOperation.isFASTQSeq(seqFile1) && SeqOperation.isFASTQSeq(seqFile2))
						   		   || (SeqOperation.isFASTASeq(seqFile1) && SeqOperation.isFASTASeq(seqFile2))
						    ){
						      System.out.println("========Filtering sequences (file No."+(i+1)+")........");
					    	  SeqFilter sf=new SeqFilter();
					    	  if(doSplitSeq) sf.setSplitStep(splitStep);
						      sf.filter(seqFile1,seqFile2,doDeDup, minSeqLen, outDir);
						      sf=null;
						  }
					    }else {
					   	  System.out.println("Warnning: No."+(i+1)+" illegal seq file, please check it!"); 
					    }
					  }
					}
					seqPairFileList=null;
					seqFileList1=null;
					seqFileList2=null;
				  }
			   }
	    	   
	    	   seqRC=null;
		   }
		   	
		   if(tmpFiles!=null){
			 for(String tmpFile: tmpFiles){
			   FileOperation.delFile(tmpFile);
			 }
		   }
		   delTmpDir(tmpDir);	
    }

}
