package org.geatools.data.structure;
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
public class  SeqCompo implements Cloneable{
	
	public boolean isActive=true;
	
	public String name="";
	public String barName=null;
	public String barSeq=null;	
	public String primerSeq=null;	
	public String primerContSeq=null;	
	public String baitTerritorySeq=null;
	public String baitTerritoryArmSeq=null;	

	/*
	public String baitSiteName=null;	
	public String baitSiteChr=null;
	public String baitSiteChrStart=null;
	public String baitSiteChrEnd=null;
	public String baitSiteStrand="+";
	*/
	public String rPrimerSeq=null;	
	public String rPrimerContSeq=null;
	public String freqCutterSeq=null;
	
	  @Override
	public Object clone() throws CloneNotSupportedException {
		SeqCompo cloned = (SeqCompo)super.clone();
		
		return cloned;
	}
	    
}
