package org.geatools;
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

import java.io.File;
import java.io.UnsupportedEncodingException;
import java.net.URLDecoder;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Calendar;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.geatools.operation.FileOperation;
import org.geatools.operation.SeqOperation;
import org.geatools.call.CallFq2Fa;
import org.geatools.call.CallQuerySeqCount;
import org.geatools.call.CallSeqFilter;
import org.geatools.call.CallSeqMut;
import org.geatools.call.CallSeqOffTarget;
import org.geatools.call.CallSeqQCFilter;
import org.geatools.call.CallSeqRecognition;
import org.geatools.call.CallSeqRetrieve;
import org.geatools.call.CallUtility;
import org.geatools.seqprocess.SeqQCFilter;
import org.geatools.seqprocess.SeqRocketConsole;

public class GEAT{
	 
	 protected static String fileSeparator = System.getProperties().getProperty("file.separator");
	 protected static String homeDir=getClassPath();
	 protected static String dataDir=homeDir+fileSeparator+"data";
	 protected static String workingDir=homeDir+fileSeparator+"working";	
	 protected static String tmpDir=homeDir+fileSeparator+"tmp";	
	 protected static List<String> tmpFiles=new ArrayList<String>();
		   
	 // basic parameters................
	 protected static String taskName=null;
	 protected static String fastq = null;
	 protected static String fasta = null;
	 protected static String fastq2 = null;
	 protected static String fasta2 = null;
	 protected static List<String> fastqList;	
	 protected static List<String> fastaList;
	 protected static List<String> fastqList2;
	 protected static List<String> fastaList2;
	 protected static String seqType=SeqOperation.SEQTYPE_SINGLEEND;
	 protected static String expSeqInfo = null;
	 protected static String expSeqInfo2 = null;	
	 protected static String outDir = null;
	 protected static String outTag=null;
	 protected static String outName=null;
	 protected static boolean isFastqOK=false;
	 protected static boolean isFastaOK=false;
	 protected static boolean isFastq2OK=false;
	 protected static boolean isFasta2OK=false;
	 protected static boolean isLibExpSeqInfoOK=false;
	 protected static boolean doSeqQCFilter=false;		
	 protected static boolean doSeqExtraction=false;
	 protected static boolean doSeqExclusion=false;	 
	 protected static SeqQCFilter seqQC=null;
	 protected static String[] seqQCOpts=null;
	 protected static SeqRocketConsole seqRC=null;
	 protected static boolean isSeqRocketOK=false;
	 protected static boolean isSeqPairRocketOK=false;	
	
	 
  public GEAT(){
	 
  }
	 
  public static Map<String, List<String>> getCommandLineParams(String[] args){
		 
		 final Map<String, List<String>> params = new HashMap<String, List<String>>();

		 List<String> options = null;
		 for (int i = 0; i < args.length; i++) {
		     final String a = args[i];

		     if (a.charAt(0) == '-') {
		         if (a.length() < 2) {
		             System.err.println("Error at argument " + a);
		             return null;
		         }

		         options = new ArrayList<String>();
		         params.put(a.substring(1), options);
		     }else if (options != null) {
		         options.add(a);
		     }else {
		         System.err.println("Illegal parameter usage");
		         return null;
		     }
		 }
		 
		 return params;
  }
	   
  public static String getClassPath() {
      String path = GEAT.class.getProtectionDomain().getCodeSource().getLocation().getPath();
      String decodedPath = path;
      try {
          decodedPath = URLDecoder.decode(path, "UTF-8");
      } catch (UnsupportedEncodingException e) {
          e.printStackTrace();
          return null;
      }

      String absolutePath = decodedPath.substring(0, decodedPath.lastIndexOf(File.separator));
      return absolutePath;
  }
  
  public static void setHomeDir(String str){
	 File dir=new File(str);
	 if(!dir.exists()) FileOperation.newFolder(str);
	 dir=null;
	 homeDir=str;	 
  }
  public static String getHomeDir(){
	 return homeDir;
  }
  
  public static void setfileSeparator(String str){
	 fileSeparator=str;
  }
  public static String getfileSeparator(){
	 return fileSeparator;
  }

  public static void setWorkingDir(String str){
	 File dir=new File(str);
	 if(!dir.exists()) FileOperation.newFolder(str);
	 dir=null;
	 workingDir=str;
  }
  public static String getWorkingDir(){
	 return workingDir;
  }
  
  public static void setDataDir(String str){
	 File dir=new File(str);
	 if(!dir.exists()) FileOperation.newFolder(str);
	 dir=null;
	 dataDir=str;
  }
  public static String getDataDir(){
	 return dataDir;
  }

  public static void setTmpDir(String str){
	 File dir=new File(str);
	 if(!dir.exists()) FileOperation.newFolder(str);
	 dir=null;
	 tmpDir=str;
  }
  public static String getTmpDir(){
	 return tmpDir;
  }

  public static void delTmpDir(String dir){
	 FileOperation.delFolder(dir);
  }
	 
  public static void main(String[] args) throws Exception {		  
	 	 
	 fileSeparator = System.getProperties().getProperty("file.separator");
	 homeDir=getClassPath();
	 dataDir=homeDir+fileSeparator+"data";
	 //String workDir = System.getProperties().getProperty("user.dir");
	 workingDir=homeDir+fileSeparator+"working";
	 String timeStamp=new SimpleDateFormat("yyyyMMdd_HHmmss").format(Calendar.getInstance().getTime());
	 String tmpDir0=homeDir+fileSeparator+"tmp";
	 tmpDir=tmpDir0+fileSeparator+timeStamp;
	
	 //System.out.println(homeDir);
	 File dir=new File(dataDir);
	 if(!dir.exists()) FileOperation.newFolder(dataDir);	
	 
	 dir=new File(workingDir);
	 if(!dir.exists()) FileOperation.newFolder(workingDir);	 

	 dir=new File(tmpDir0);
	 if(!dir.exists()) FileOperation.newFolder(tmpDir0);
	 dir=new File(tmpDir);
	 if(!dir.exists()) FileOperation.newFolder(tmpDir);
	 dir=null;
	   
	 Map<String, List<String>> params=getCommandLineParams(args);
	 
	 taskName="Utility";
	 
	 if(params.get("task")!=null){		
		if(params.get("task").size()>0){
		  taskName=params.get("task").get(0).trim();		 
		}
	 }
	 
   
	 if(taskName.equalsIgnoreCase("SeqRecognition")){
		 CallSeqRecognition.doWork(args); 
     }else if(taskName.equalsIgnoreCase("SeqMut")){
		 CallSeqMut.doWork(args);
	 }else if(taskName.equalsIgnoreCase("SeqOffTarget")){
	     CallSeqOffTarget.doWork(args);	  
	 }else if(taskName.equalsIgnoreCase("SeqQCFilter")){
		 CallSeqQCFilter.doWork(args);			    	
	 }else if(taskName.equalsIgnoreCase("Fq2Fa")){
		 CallFq2Fa.doWork(args);			    	
	 }else if(taskName.equalsIgnoreCase("QuerySeqCount")){
		 CallQuerySeqCount.doWork(args);
     }else if(taskName.equalsIgnoreCase("SeqFilter")){
		 CallSeqFilter.doWork(args);			    	
	 }else if(taskName.equalsIgnoreCase("SeqRetrieve")){
		 CallSeqRetrieve.doWork(args);			    	
	 }else if(taskName.equalsIgnoreCase("Utility")){
		 CallUtility.doWork(args);			 	  
	 }else{// if -task
		 System.err.println("Error: '-task' is invalid");				
	 }
	  
	 delTmpDir(tmpDir);
  }
}