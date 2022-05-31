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
import org.geatools.data.load.LoadChrInfo;
import org.geatools.data.structure.ChrInfo;
import org.geatools.data.structure.SeqAlignSite;
import org.geatools.data.structure.SeqInfo;
import org.geatools.operation.FileOperation;
import org.geatools.operation.SeqOperation;
import org.geatools.seqprocess.SeqOffTarget;

public class CallSeqOffTarget extends GEAT{
	
	public static void doWork( String[] args){
		 
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

		 String targetSeqFile = null;
		 String targetSeqStr;
		 String refGenome = "mm9";
		 float seqMinIdentity = 80.0f;
		 float seqMinAlignLenRatio = 0.7f; 
		 float seqMaxMismatchRatio=0.25f;
		 int maxGapNum=2;
		 String outDir = null;	 
	     int siteDownLen=6;
	     int siteUpLen=0;
	     
		 boolean doIt=false;
		 String taskName;
		 //boolean isTargetSeqFileOK=false;		 

		 Map<String, List<String>> params=getCommandLineParams(args);
		 if(params.get("task")!=null){		
			doIt=false;	
			if(params.get("task").size()>0){
			  taskName=params.get("task").get(0).trim();
			  if(taskName!=null){
			    if(taskName.equalsIgnoreCase("SeqOffTarget")){
			   	  doIt=true;	    	    
			    }else{
			      doIt=false;
			      System.err.println("Error: '-task' is invalid");
				  return;	
			    }
			  }
			}
		 }else{
		   System.err.println("Error: '-task' is invalid");
		   return;		
		 }
		 
		 if(doIt){
		     //####### set target Seq #######
		     if(params.get("targetSeq")!=null){
				  if(params.get("targetSeq").size()>0){
					 targetSeqStr=params.get("targetSeq").get(0).trim();
				     if(targetSeqStr==null){
				    	targetSeqFile=tmpDir+fileSeparator+"YourTargetSeq.fna";
				    	SeqOperation.saveSeqAsFASTA(targetSeqStr,"YourTargetSeq",targetSeqFile);
				     }else{
				    	System.err.println("Illegal '-targetSeq' parameter usage :(");
						return;
				     }
				  }else{
					 System.err.println("Illegal '-targetSeq' parameter usage :(");
					 return;
				  }
		     }	
			
			//####### set target Seq fasta #######
			 if(params.get("targetSeqFASTA")!=null){		
			   if(params.get("targetSeqFASTA").size()>0){
				 targetSeqFile=params.get("targetSeqFASTA").get(0).trim();
			     if(targetSeqFile!=null){
			       File f = new File(targetSeqFile);
				   if(!f.exists()){
					 System.err.println("'-targetSeqFASTA' file you provided doesn't exist :(");
				     return;
				   }
				   f=null;
				 }
			   }else{
				 System.err.println("You didn't provide '-targetSeqFASTA' file :(");
				 return;
			   }
			 }	
		     
			 //####### set refGenome assembly #######
			 if(params.get("refGenome")!=null){		  	
			   if(params.get("refGenome").size()>0){
			      refGenome=params.get("refGenome").get(0).trim();
				  if(refGenome==null){			   	
					System.err.println("Illegal '-refGenome' parameter usage :(");
				    return;
				  }
			   }else{
				  System.err.println("Illegal '-refGenome' parameter usage :(");
				  return;
			   }
			 }  
	         
			 //####### set seq minmum align length to target seq length ratio #######
			 if(params.get("minAlignLenRatio")!=null){		  	
			   if(params.get("minAlignLenRatio").size()>0){
				  String tmp=params.get("minAlignLenRatio").get(0).trim();
				  if(tmp!=null){
					try{
					   seqMinAlignLenRatio=Float.parseFloat(tmp);
					   if(seqMinAlignLenRatio>1 || seqMinAlignLenRatio<0){ 
						 System.err.println("Illegal '-minAlignLenRatio' parameter usage :(");
						 return;
					   }	
					}catch(Exception e){
					   System.err.println("Illegal '-minAlignLenRatio' parameter usage :(");
					   return;
				    } 
				  }else{
					System.err.println("Illegal '-minAlignLenRatio' parameter usage :(");
				    return;
				  }
			   }else{
				  System.err.println("Illegal '-minAlignLenRatio' parameter usage :(");
				  return;
			   }
			 }
			 
			 //####### set seq minmum align identity #######
			 if(params.get("minIdentity")!=null){		  	
			   if(params.get("minIdentity").size()>0){
				  String tmp=params.get("minIdentity").get(0).trim();
				  if(tmp!=null){
					try{
					 seqMinIdentity=Float.parseFloat(tmp);
					 if(seqMinIdentity>=0 && seqMinIdentity<=1){ 
						 seqMinIdentity=seqMinIdentity*100.0f;
					 }else if(seqMinIdentity>100 || seqMinIdentity<0){ 
						 System.err.println("Illegal '-minIdentity' parameter usage :(");
						 return;
					 }		 
					}catch(Exception e){
					   System.err.println("Illegal '-minIdentity' parameter usage :(");
					   return;
				    } 
				  }else{
					System.err.println("Illegal '-minIdentity' parameter usage :(");
				    return;
				  }
			   }else{
				  System.err.println("Illegal '-minIdentity' parameter usage :(");
				  return;
			   }
			 }
			 
			 //####### set seq max mismatch ratio #######
			 if(params.get("maxMismatchRatio")!=null){		  	
			   if(params.get("maxMismatchRatio").size()>0){
				  String tmp=params.get("maxMismatchRatio").get(0).trim();
				  if(tmp!=null){
					try{
					  seqMaxMismatchRatio=Float.parseFloat(tmp);
					  if(seqMaxMismatchRatio>1 || seqMaxMismatchRatio<0){ 
						System.err.println("Illegal '-maxMismatchRatio' parameter usage :(");
						return;
					  }	
					}catch(Exception e){
					   System.err.println("Illegal '-maxMismatchRatio' parameter usage :(");
					   return;
				    } 
				  }else{
					System.err.println("Illegal '-maxMismatchRatio' parameter usage :(");
				    return;
				  }
			   }else{
				  System.err.println("Illegal '-maxMismatchRatio' parameter usage :(");
				  return;
			   }
			 }
			 
			 //####### set seq max gap number #######
			 if(params.get("maxGapNum")!=null){		  	
			   if(params.get("maxGapNum").size()>0){
				  String tmp=params.get("maxGapNum").get(0).trim();
				  if(tmp!=null){
					try{
					  maxGapNum=Integer.parseInt(tmp);
					  if(maxGapNum<0){ 
						maxGapNum=0;
					  }	
					}catch(Exception e){
					   System.err.println("Illegal '-maxGapNum' parameter usage :(");
					   return;
				    } 
				  }else{
					System.err.println("Illegal '-maxGapNum' parameter usage :(");
				    return;
				  }
			   }else{
				  System.err.println("Illegal '-maxGapNum' parameter usage :(");
				  return;
			   }
			 }  
		 	 
			//####### set output #######
			 boolean  doOutput=false;
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
			   }
			 }
			 if(!doOutput){			   
				 outDir=workingDir;
				 dir=new File(outDir);
				 if(!dir.exists()) FileOperation.newFolder(outDir);	
		         dir=null;
				 doOutput=true;   	
			 }	
		     
		 	 System.out.println("Loading Chr......");	
		 	 String genomePath="data"+fileSeparator+"genome_info"+fileSeparator+refGenome+fileSeparator;
		 	 String chrInfoFile=genomePath+fileSeparator+refGenome+"_chromInfo.txt";
		 	 String chrSeqPath=genomePath+fileSeparator+refGenome+"_chromfa"+fileSeparator;
		 	 List<ChrInfo> chrInfoList;
		     chrInfoList=LoadChrInfo.getChrInfo(chrInfoFile);
		     LoadChrInfo.sortChrByNum(chrInfoList);    
		      	 	
		 	 String chrName=""; 	
		 	 String chrFile="";
		 	 for(int i=0;i<chrInfoList.size(); i++){
		 	   chrName=chrInfoList.get(i).name;
		 	   chrFile=chrSeqPath+fileSeparator+chrName+".fa";
		 	   chrInfoList.get(i).chrFaFile=chrFile;
		 	 }
		 	 
		 	 System.out.println("Scanning Off-target Site......");	
		 	 SeqOffTarget OTS=new SeqOffTarget();
		 	 OTS.setTmpDir(tmpDir);
		     List<SeqInfo> targetSeqList=SeqOperation.getFASTASeqObj(targetSeqFile);
		     int queryCount=0;
		 	 for(SeqInfo targetSeq:targetSeqList){
		 		 queryCount++;
			 	 List<SeqAlignSite> chrOfftargetSites;
			 	 List<ArrayList<String>> offTargetList=new ArrayList<ArrayList<String>>();
			 	 ArrayList<String> offTarget=new ArrayList<String>();
			 	 if(targetSeq.seqName==null || targetSeq.seqName==""){
			 		targetSeq.seqName="TargetSeq"+queryCount;
			 	 }
			 	 System.out.println("### For target sequence "+targetSeq.seqName+"......");	
			 	 String targetSeqFaFile=tmpDir+targetSeq.seqName+".fna";
			 	 SeqOperation.saveSeqObjAsFASTA(targetSeq, targetSeqFaFile);
			 	 String genomeChrSeqFile="";
			 	 
			  	 for(int i=0;i<chrInfoList.size(); i++){
			  	   System.out.println("Scanning "+chrInfoList.get(i).name+"......");	
			  	   genomeChrSeqFile=chrInfoList.get(i).chrFaFile;
			 	   
			 	   chrOfftargetSites=OTS.getOffTargetSite(targetSeqFaFile,genomeChrSeqFile,
			 			   seqMinIdentity,seqMinAlignLenRatio,seqMaxMismatchRatio,maxGapNum); 
			 	   
			 	   //sortSiteByScore(chrOfftargetSites,true);
			 	   OTS.sortSiteByMismatchNum(chrOfftargetSites,false);
			 	  
			 	   List<String> chrLineSeq=SeqOperation.getGenomeChrLineSeq(genomeChrSeqFile);
			 	   
			 	   SeqOperation.setAlignSiteSeq(chrOfftargetSites,chrLineSeq,siteDownLen,siteUpLen);
			 	   List<Integer> targetSiteIdx=new ArrayList<Integer>();
			 	   for(int j=0;j<chrOfftargetSites.size();j++){	
			 		  if(chrOfftargetSites.get(j).name.equalsIgnoreCase("OffTarget")){
				 		 offTarget=new ArrayList<String>();	 
				 		 offTarget.add(chrOfftargetSites.get(j).chr); 	 
				 		 offTarget.add(Integer.toString(chrOfftargetSites.get(j).chrStart));
				 		 offTarget.add(Integer.toString(chrOfftargetSites.get(j).chrEnd));
				 		 offTarget.add(chrOfftargetSites.get(j).strand);
				 		 offTarget.add(chrOfftargetSites.get(j).name);
				 		 offTarget.add(Double.toString(chrOfftargetSites.get(j).score));
				 		 offTarget.add(Integer.toString(chrOfftargetSites.get(j).mismatchNum));
				 		 offTarget.add(chrOfftargetSites.get(j).seq);
				 		 offTargetList.add(offTarget);
				 	     offTarget=null;
			 		  }else if(chrOfftargetSites.get(j).name.equalsIgnoreCase("Target")){
			 			 targetSiteIdx.add(j);
			 		  }
			 	   }
			 	   
			 	   for(int t=0;t<targetSiteIdx.size();t++){	
			 	         int j=targetSiteIdx.get(t);
				 		 offTarget=new ArrayList<String>();	 
				 		 offTarget.add(chrOfftargetSites.get(j).chr); 	 
				 		 offTarget.add(Integer.toString(chrOfftargetSites.get(j).chrStart));
				 		 offTarget.add(Integer.toString(chrOfftargetSites.get(j).chrEnd));
				 		 offTarget.add(chrOfftargetSites.get(j).strand);
				 		 offTarget.add(chrOfftargetSites.get(j).name);
				 		 offTarget.add(Double.toString(chrOfftargetSites.get(j).score));
				 		 offTarget.add(Integer.toString(chrOfftargetSites.get(j).mismatchNum));
				 		 offTarget.add(chrOfftargetSites.get(j).seq);
				 		 offTargetList.add(t,offTarget);
				 	     offTarget=null;			 		  
			 	   }
			 	   targetSiteIdx=null;
			 	}
			 	
			 	String outFile=outDir+fileSeparator+refGenome
			 			+targetSeq.seqName+".OffTargetSite_IP"+seqMinIdentity
			 			+"LR"+seqMinAlignLenRatio+"MR"+seqMaxMismatchRatio
			 			+"G"+maxGapNum+".bed";
			 	FileOperation.saveMatrixList(offTargetList, outFile);
			
			 	offTargetList=null;
		 	}
		 	OTS.delTmpDir(tmpDir);
		 	OTS=null;
		 	targetSeqList=null;
		 	delTmpDir(tmpDir);
	    }// if doIt  
	  }
}
