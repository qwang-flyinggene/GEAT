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
import java.util.List;
import java.util.Map;

import org.geatools.GEAT;
import org.geatools.NumberCheck;
import org.geatools.data.load.LoadChrInfo;
import org.geatools.data.load.LoadChrSite;
import org.geatools.data.structure.ChrInfo;
import org.geatools.operation.FileOperation;
import org.geatools.operation.SeqOperation;


public class CallSeqRetrieve extends GEAT{
	
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

		 String genomeRef = null;
		 String chrSiteFile="";
		 String fileOfChrSiteFile="";
		 String siteChr="";
		 int siteStart=0;
		 int siteEnd=0;
		 String siteStrand="+";
	     int siteDownLen=0;
	     int siteUpLen=0;
	     boolean doRCForMinusStrand=true;
	     boolean saveChrSep=false;
	     boolean saveChrTogether=true;
	     
		 boolean doIt=false;
		 String taskName;
		 //boolean isTargetSeqFileOK=false;
		 			
		 Map<String, List<String>> params=getCommandLineParams(args);
		 if(params.get("task")!=null){		
			doIt=false;	
			if(params.get("task").size()>0){
			  taskName=params.get("task").get(0).trim();
			  if(taskName!=null){
			    if(taskName.equalsIgnoreCase("SeqRetrieve")){
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

			 
			 //####### set refGenome assembly #######
			 if(params.get("ref")!=null){		  	
			   if(params.get("ref").size()>0){
				  genomeRef=params.get("ref").get(0).trim();
				  if(genomeRef==null){			   	
					System.err.println("Illegal '-ref' parameter usage :(");
				    return;
				  }
			   }else{
				  System.err.println("Illegal '-ref' parameter usage :(");
				  return;
			   }
			 }
			 
			 //####### set refGenome assembly #######
			 if(params.get("RCMinus")!=null){	
			   doRCForMinusStrand=true;
			   if(params.get("RCMinus").size()>0){
			     String str=params.get("RCMinus").get(0).trim();
				 if(str.equalsIgnoreCase("no") 
					  || str.equalsIgnoreCase("n")
					  || str.equalsIgnoreCase("false")
					  || str.equalsIgnoreCase("f")){			   	
					 
					 doRCForMinusStrand=false;
				 }
			   }
			 }
			 
			 if(params.get("saveChrSep")!=null){	
				saveChrSep=true;
				if(params.get("saveChrSep").size()>0){
				   String str=params.get("saveChrSep").get(0).trim();
				   if(str.equalsIgnoreCase("no") 
					  || str.equalsIgnoreCase("n")
					  || str.equalsIgnoreCase("false")
					  || str.equalsIgnoreCase("f")){			   	
						 
					   saveChrSep=false;
				   }
				}
			 }
			 
			 if(params.get("saveChrTogether")!=null){	
				saveChrTogether=true;
				if(params.get("saveChrTogether").size()>0){
				   String str=params.get("saveChrTogether").get(0).trim();
				   if(str.equalsIgnoreCase("no") 
					  || str.equalsIgnoreCase("n")
					  || str.equalsIgnoreCase("false")
					  || str.equalsIgnoreCase("f")){			   	
						 
					   saveChrTogether=false;
				   }
				}
			 }
			 
			 boolean isChrSiteOK=false;
			 //####### set upStream #######
			 if(params.get("chrSite")!=null){		  	
			   if(params.get("chrSite").size()>=3){
				  isChrSiteOK=true;
				  String str=params.get("chrSite").get(0).trim();
				  if(str.matches(LoadChrInfo.getChrRegex())) {
					siteChr=str;
				  }else{			   	
					System.err.println("Illegal '-chrSite chr start end' parameter usage, 'chr' name is incorrect!");
					isChrSiteOK=false;
					//return;
				  }
				  
				  str=params.get("chrSite").get(1).trim();
				  if(NumberCheck.isPositiveInteger(str)) {
					siteStart=Integer.parseInt(str);
				  }else{			   	
					System.err.println("Illegal '-chrSite chr start end' parameter usage, 'start' is incorrect!");
					isChrSiteOK=false;
					//return;
				  }
				  
				  str=params.get("chrSite").get(2).trim();
				  if(NumberCheck.isPositiveInteger(str)) {
				    siteEnd=Integer.parseInt(str);
				  }else{			   	
					System.err.println("Illegal '-chrSite chr start end' parameter usage, 'end' is incorrect!");
					isChrSiteOK=false;
					//return;
				  }			  
				 				  
				  str=params.get("chrSite").get(3).trim();
				  if(str.equals("plus")) {
				    siteStrand="+";
				  }else if(str.equals("minus")){
					siteStrand="-"; 
				  }			  
				  
			   }else{
				  System.err.println("Illegal '-chrSite' parameter usage, it should be inputted by '-chrSite chr start end strand(optional)!");
				  //return;
			   }
			 }  
			 
			 boolean isChrSiteListOK=false;
			 //####### set chrSiteList #######
			 if(params.get("chrSitesFile")!=null){		  	
			   if(params.get("chrSitesFile").size()>0){
			      chrSiteFile=params.get("chrSitesFile").get(0).trim();
			      isChrSiteListOK=true;
				  if(chrSiteFile==null || !new File(chrSiteFile).exists()){
					isChrSiteListOK=false;
					System.err.println("Illegal '-chrSitesFile' parameter usage or provided file does not exist!");					
					//return;
				  }				 
			   }else{
				  System.err.println("Illegal '-chrSitesFile' parameter usage!");
				  //return;
			   }
			 } 
			 
			 boolean isFilesOfChrSiteListOK=false;
			 //####### set chrSiteList #######
			 if(params.get("list_chrSitesFile")!=null){		  	
			   if(params.get("list_chrSitesFile").size()>0){
				  fileOfChrSiteFile=params.get("list_chrSitesFile").get(0).trim();
				  isFilesOfChrSiteListOK=true;
				  if(fileOfChrSiteFile==null || !new File(fileOfChrSiteFile).exists()){
					isFilesOfChrSiteListOK=false;
					System.err.println("Illegal '-list_chrSitesFile' parameter usage or provided file does not exist!");					
					//return;
				  }				 
			   }else{
				  System.err.println("Illegal '-list_chrSitesFile' parameter usage!");
				  //return;
			   }
			 } 

			 
			 if(!isChrSiteOK && !isChrSiteListOK && !isFilesOfChrSiteListOK) {
				 System.err.println("You must correctly set at least one parameter of '-chrSiteList', '-chrSite' or 'multi_chrSiteList' :(");
				 return;
			 }

			 //####### set upStream #######
			 if(params.get("upStream")!=null){		  	
			   if(params.get("upStream").size()>0){
				  String str=params.get("upStream").get(0).trim();
				  if(NumberCheck.isPositiveInteger(str)) {
					siteUpLen=Integer.parseInt(str);
				  }else{			   	
					System.err.println("Illegal '-upStream' parameter usage!");
				    return;
				  }
			   }else{
				  System.err.println("Illegal '-upStream' parameter usage!");
				  return;
			   }
			 }  
			 
			 //####### set downStream #######
			 if(params.get("downStream")!=null){		  	
			   if(params.get("downStream").size()>0){
				  String str=params.get("downStream").get(0).trim();
				  if(NumberCheck.isPositiveInteger(str)) {
					siteDownLen=Integer.parseInt(str);
				  }else{			   	
					System.err.println("Illegal '-downStream' parameter usage!");
				    return;
				  }
			   }else{
				  System.err.println("Illegal '-downStream' parameter usage!");
				  return;
			   }
			 }  
		     
		 	 //System.out.println("Loading Chromosomes......");	
			 if(genomeRef==null) {
				System.err.println("Error! You did not set '-ref' parameter, please set it and run again!");
				return;
			 }
		 	 String genomePath="data"+fileSeparator+"genome_info"+fileSeparator+genomeRef+fileSeparator;
		 	 String chrInfoFile=genomePath+fileSeparator+genomeRef+"_chromInfo.txt";
		 	 String chrSeqPath=genomePath+fileSeparator+genomeRef+"_chromfa"+fileSeparator;
		 	 
		 	 if(!new File(chrInfoFile).exists()) {
				System.err.println("Error! The 'chromInfo' file ["+chrInfoFile+"] does not exist, please check it!");
				return; 
		 	 }
		 	 if(!new File(chrSeqPath).exists()) {
				System.err.println("Error! The 'chromfa' directories ["+chrSeqPath+"] does not exist, please check it!");
				return; 
		 	 }

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
		 	 
			 String genomeChrSeqFile="";
			 List<String> chrLineSeq;	
		 	 if(isChrSiteOK) {		 		 
			   for(int i=0;i<chrInfoList.size(); i++){
				 if(siteChr.equalsIgnoreCase(chrInfoList.get(i).name)) {
				   System.out.println("Scanning "+chrInfoList.get(i).name+"......");	
				   genomeChrSeqFile=chrInfoList.get(i).chrFaFile;		 	   
        		   chrLineSeq=SeqOperation.getGenomeChrLineSeq(genomeChrSeqFile);
				   String seq=SeqOperation.getChrSiteSeq(chrLineSeq,siteStart,siteEnd,siteStrand,siteDownLen,siteUpLen);					 	   
				   System.out.println("Retrieved site: strand("+siteStrand+") "
				    		+ siteChr+":"+siteStart+"-"+siteEnd
				    		+" plus "+siteUpLen+"bp upstream and "+siteDownLen+"bp downstream");
				   if(siteStrand.equals("-")) seq=SeqOperation.reverseComplement(seq);
				   System.out.println("Sequence 5'->3' strand("+siteStrand+"): "+ seq);
				   chrLineSeq=null;
				   break;
				 }
			   }
		 
		 	 }		 	 
		 	 
		 	 if(isChrSiteListOK) {
		 	   System.out.println("Working on chr sites file ["+chrSiteFile+"]......");	
		 	   if(outTag==null) outTag="chrSiteList";
		 	   if(outName==null) { 
		 		  outName=chrSiteFile.substring(chrSiteFile.lastIndexOf(fileSeparator)+1,chrSiteFile.lastIndexOf("."))
		 		       +".RetrievedSeq.fna";
		 	   }
		 	   if(outDir==null) outDir=chrSiteFile.substring(0,chrSiteFile.lastIndexOf(fileSeparator));
		 	   
		 	   String outFile=outDir+fileSeparator+outName;
		 	   String outChrFile;
		 	   LoadChrSite.setChrSites(chrSiteFile,chrInfoList);		 	 
			   for(ChrInfo chrInfo:chrInfoList){
			     System.out.println("Scanning "+chrInfo.name+"......");	
			     if(chrInfo.siteList!=null && chrInfo.siteList.size()>0){
		   	       genomeChrSeqFile=chrInfo.chrFaFile;	
			 	   chrLineSeq=SeqOperation.getGenomeChrLineSeq(genomeChrSeqFile);			 	   
			 	   SeqOperation.setChrSiteSeq(chrInfo.siteList,chrLineSeq,siteDownLen,siteUpLen,doRCForMinusStrand);
			 	   if(saveChrSep) {
			 	     outChrFile=outFile.substring(0,outFile.lastIndexOf("."))+"."+chrInfo.name+".fna";
			 	     LoadChrSite.saveChrSiteSeqAsFASTA(chrInfo.siteList, outChrFile);
			 	   }
			 	   chrLineSeq=null;  
			     }
			   }
			   
			   if(saveChrTogether) LoadChrSite.saveChrSiteSeqAsFASTA(LoadChrSite.getChrSites(chrInfoList),outFile);
			   System.out.println("Sequence retrive done!");
		 	 }
		 	 
		 	 if(isFilesOfChrSiteListOK) {
		 		 List<String> fileListOfChrSiteFile=FileOperation.getRowsOfFile(fileOfChrSiteFile);
		 		 int i=0;
		 		 for(String siteFile: fileListOfChrSiteFile) {	
		 		   i++;
		 		   System.out.println("Working on No. "+i+" chr sites file ["+siteFile+"]......");	
			 	  
			 	   outName=siteFile.substring(siteFile.lastIndexOf(fileSeparator)+1,siteFile.lastIndexOf("."))
			 		       +".RetrievedSeq.fna";
			 	   
			 	   if(outDir==null) outDir=siteFile.substring(0,siteFile.lastIndexOf(fileSeparator));
			 	   
			 	   String outFile=outDir+fileSeparator+outName;
			 	   String outChrFile;
			 	   LoadChrSite.setChrSites(siteFile,chrInfoList);		 	 
				   for(ChrInfo chrInfo:chrInfoList){
				     System.out.println("Scanning "+chrInfo.name+"......");	
				     if(chrInfo.siteList!=null && chrInfo.siteList.size()>0){
			   	       genomeChrSeqFile=chrInfo.chrFaFile;	
				 	   chrLineSeq=SeqOperation.getGenomeChrLineSeq(genomeChrSeqFile);			 	   
				 	   SeqOperation.setChrSiteSeq(chrInfo.siteList,chrLineSeq,siteDownLen,siteUpLen,doRCForMinusStrand);
				 	   if(saveChrSep) {
				 	     outChrFile=outFile.substring(0,outFile.lastIndexOf("."))+"."+chrInfo.name+".fna";
				 	     LoadChrSite.saveChrSiteSeqAsFASTA(chrInfo.siteList, outChrFile);
				 	   }
				 	   chrLineSeq=null;  
				     }
				   }
				   
				   if(saveChrTogether) LoadChrSite.saveChrSiteSeqAsFASTA(LoadChrSite.getChrSites(chrInfoList),outFile);
		 		 }
		 		 fileListOfChrSiteFile=null;
		 		 System.out.println("Sequence retrive done!");
		 	 }		 	

		 	 delTmpDir(tmpDir);
	    }// if doIt  
	  }
}
