/*
 * Copyright 1997, Regents of the University of Minnesota
 *
 * rename.h
 *
 * This file contains header files
 *
 * Started 10/2/97
 * George
 *
 * $Id: rename.h,v 1.2 2004/01/31 00:38:31 bates Exp $
 *
 */

/* balance.c */
#define Balance2Way			__Balance2Way
#define Bnd2WayBalance			__Bnd2WayBalance
#define General2WayBalance		__General2WayBalance


/* bucketsort.c */
#define BucketSortKeysInc		__BucketSortKeysInc


/* ccgraph.c */
#define CreateCoarseGraph		__CreateCoarseGraph
#define CreateCoarseGraphNoMask		__CreateCoarseGraphNoMask
#define CreateCoarseGraph_NVW 		__CreateCoarseGraph_NVW
#define SetUpCoarseGraph		__SetUpCoarseGraph
#define ReAdjustMemory			__ReAdjustMemory


/* coarsen.c */
#define Coarsen2Way			__Coarsen2Way


/* compress.c */
#define CompressGraph			__CompressGraph
#define PruneGraph			__PruneGraph


/* debug.c */
#define ComputeCut			__ComputeCut
#define CheckBnd			__CheckBnd
#define CheckBnd2			__CheckBnd2
#define CheckNodeBnd			__CheckNodeBnd
#define CheckRInfo			__CheckRInfo
#define CheckNodePartitionParams	__CheckNodePartitionParams
#define IsSeparable			__IsSeparable


/* estmem.c */
#define EstimateCFraction		__EstimateCFraction
#define ComputeCoarseGraphSize		__ComputeCoarseGraphSize


/* fm.c */
#define FM_2WayEdgeRefine		__FM_2WayEdgeRefine


/* fortran.c */
#define Change2CNumbering		__Change2CNumbering
#define Change2FNumbering		__Change2FNumbering
#define Change2FNumbering2		__Change2FNumbering2
#define Change2FNumberingOrder		__Change2FNumberingOrder
#define ChangeMesh2CNumbering		__ChangeMesh2CNumbering
#define ChangeMesh2FNumbering		__ChangeMesh2FNumbering
#define ChangeMesh2FNumbering2		__ChangeMesh2FNumbering2


/* graph.c */
#define SetUpGraph			__SetUpGraph
#define SetUpGraphKway 			__SetUpGraphKway
#define SetUpGraph2			__SetUpGraph2
#define VolSetUpGraph			__VolSetUpGraph
#define RandomizeGraph			__RandomizeGraph
#define IsConnectedSubdomain		__IsConnectedSubdomain
#define IsConnected			__IsConnected
#define IsConnected2			__IsConnected2
#define FindComponents			__FindComponents


/* initpart.c */
#define Init2WayPartition		__Init2WayPartition
#define InitSeparator			__InitSeparator
#define GrowBisection			__GrowBisection
#define GrowBisectionNode		__GrowBisectionNode
#define RandomBisection			__RandomBisection


/* kmetis.c */
#define MlevelKWayPartitioning		__MlevelKWayPartitioning


/* kvmetis.c */
#define MlevelVolKWayPartitioning	__MlevelVolKWayPartitioning


/* kwayfm.c */
#define Random_KWayEdgeRefine		__Random_KWayEdgeRefine
#define Greedy_KWayEdgeRefine		__Greedy_KWayEdgeRefine
#define Greedy_KWayEdgeBalance		__Greedy_KWayEdgeBalance


/* kwayrefine.c */
#define RefineKWay			__RefineKWay
#define AllocateKWayPartitionMemory	__AllocateKWayPartitionMemory
#define ComputeKWayPartitionParams	__ComputeKWayPartitionParams
#define ProjectKWayPartition		__ProjectKWayPartition
#define IsBalanced			__IsBalanced
#define ComputeKWayBoundary		__ComputeKWayBoundary
#define ComputeKWayBalanceBoundary	__ComputeKWayBalanceBoundary


/* kwayvolfm.c */
#define Random_KWayVolRefine		__Random_KWayVolRefine
#define Random_KWayVolRefineMConn	__Random_KWayVolRefineMConn
#define Greedy_KWayVolBalance		__Greedy_KWayVolBalance
#define Greedy_KWayVolBalanceMConn	__Greedy_KWayVolBalanceMConn
#define KWayVolUpdate			__KWayVolUpdate
#define ComputeKWayVolume		__ComputeKWayVolume
#define ComputeVolume			__ComputeVolume
#define CheckVolKWayPartitionParams	__CheckVolKWayPartitionParams
#define ComputeVolSubDomainGraph	__ComputeVolSubDomainGraph
#define EliminateVolSubDomainEdges	__EliminateVolSubDomainEdges


/* kwayvolrefine.c */
#define RefineVolKWay			__RefineVolKWay
#define AllocateVolKWayPartitionMemory	__AllocateVolKWayPartitionMemory
#define ComputeVolKWayPartitionParams	__ComputeVolKWayPartitionParams
#define ComputeKWayVolGains		__ComputeKWayVolGains
#define ProjectVolKWayPartition		__ProjectVolKWayPartition
#define ComputeVolKWayBoundary		__ComputeVolKWayBoundary
#define ComputeVolKWayBalanceBoundary	__ComputeVolKWayBalanceBoundary


/* match.c */
#define Match_RM			__Match_RM
#define Match_RM_NVW			__Match_RM_NVW
#define Match_HEM			__Match_HEM
#define Match_SHEM			__Match_SHEM


/* mbalance.c */
#define MocBalance2Way			__MocBalance2Way
#define MocGeneral2WayBalance		__MocGeneral2WayBalance


/* mbalance2.c */
#define MocBalance2Way2			__MocBalance2Way2
#define MocGeneral2WayBalance2		__MocGeneral2WayBalance2
#define SelectQueue3			__SelectQueue3


/* mcoarsen.c */
#define MCCoarsen2Way			__MCCoarsen2Way


/* memory.c */
#define AllocateWorkSpace		__AllocateWorkSpace
#define FreeWorkSpace			__FreeWorkSpace
#define WspaceAvail			__WspaceAvail
#define idxwspacemalloc			__idxwspacemalloc
#define idxwspacefree			__idxwspacefree
#define fwspacemalloc			__fwspacemalloc
#define CreateGraph			__CreateGraph
#define InitGraph			__InitGraph
#define FreeGraph			__FreeGraph


/* mesh.c */
#define TRIDUALMETIS			__TRIDUALMETIS
#define TETDUALMETIS			__TETDUALMETIS
#define HEXDUALMETIS			__HEXDUALMETIS
#define TRINODALMETIS			__TRINODALMETIS
#define TETNODALMETIS			Metis_TETNODALMETIS
#define HEXNODALMETIS			Metis_HEXNODALMETIS


/* mfm.c */
#define MocFM_2WayEdgeRefine		Metis_MocFM_2WayEdgeRefine
#define SelectQueue			Metis_SelectQueue
#define BetterBalance			Metis_BetterBalance
#define Compute2WayHLoadImbalance	Metis_Compute2WayHLoadImbalance
#define Compute2WayHLoadImbalanceVec	Metis_Compute2WayHLoadImbalanceVec


/* mfm2.c */
#define MocFM_2WayEdgeRefine2		Metis_MocFM_2WayEdgeRefine2
#define SelectQueue2			Metis_SelectQueue2
#define IsBetter2wayBalance		Metis_IsBetter2wayBalance


/* mincover.c */
#define MinCover			Metis_MinCover
#define MinCover_Augment		Metis_MinCover_Augment
#define MinCover_Decompose		Metis_MinCover_Decompose
#define MinCover_ColDFS			Metis_MinCover_ColDFS
#define MinCover_RowDFS			Metis_MinCover_RowDFS


/* minitpart.c */
#define MocInit2WayPartition		Metis_MocInit2WayPartition
#define MocGrowBisection		Metis_MocGrowBisection
#define MocRandomBisection		Metis_MocRandomBisection
#define MocInit2WayBalance		Metis_MocInit2WayBalance
#define SelectQueueoneWay		Metis_SelectQueueoneWay


/* minitpart2.c */
#define MocInit2WayPartition2		Metis_MocInit2WayPartition2
#define MocGrowBisection2		Metis_MocGrowBisection2
#define MocGrowBisectionNew2		Metis_MocGrowBisectionNew2
#define MocInit2WayBalance2		Metis_MocInit2WayBalance2
#define SelectQueueOneWay2		Metis_SelectQueueOneWay2


/* mkmetis.c */
#define MCMlevelKWayPartitioning	Metis_MCMlevelKWayPartitioning


/* mkwayfmh.c */
#define MCRandom_KWayEdgeRefineHorizontal	Metis_MCRandom_KWayEdgeRefineHorizontal
#define MCGreedy_KWayEdgeBalanceHorizontal	Metis_MCGreedy_KWayEdgeBalanceHorizontal
#define AreAllHVwgtsBelow			Metis_AreAllHVwgtsBelow
#define AreAllHVwgtsAbove			Metis_AreAllHVwgtsAbove
#define ComputeHKWayLoadImbalance		Metis_ComputeHKWayLoadImbalance
#define MocIsHBalanced				Metis_MocIsHBalanced
#define IsHBalanceBetterFT			Metis_IsHBalanceBetterFT
#define IsHBalanceBetterTT			Metis_IsHBalanceBetterTT


/* mkwayrefine.c */
#define MocRefineKWayHorizontal		Metis_MocRefineKWayHorizontal
#define MocAllocateKWayPartitionMemory	Metis_MocAllocateKWayPartitionMemory
#define MocComputeKWayPartitionParams	Metis_MocComputeKWayPartitionParams
#define MocProjectKWayPartition		Metis_MocProjectKWayPartition
#define MocComputeKWayBalanceBoundary	Metis_MocComputeKWayBalanceBoundary


/* mmatch.c */
#define MCMatch_RM			Metis_MCMatch_RM
#define MCMatch_HEM			Metis_MCMatch_HEM
#define MCMatch_SHEM			Metis_MCMatch_SHEM
#define MCMatch_SHEBM			Metis_MCMatch_SHEBM
#define MCMatch_SBHEM			Metis_MCMatch_SBHEM
#define BetterVBalance			Metis_BetterVBalance
#define AreAllVwgtsBelowFast		Metis_AreAllVwgtsBelowFast


/* mmd.c */
#define genmmd				Metis_genmmd
#define mmdelm				Metis_mmdelm
#define mmdint				Metis_mmdint
#define mmdnum				Metis_mmdnum
#define mmdupd				Metis_mmdupd


/* mpmetis.c */
#define MCMlevelRecursiveBisection	Metis_MCMlevelRecursiveBisection
#define MCHMlevelRecursiveBisection	Metis_MCHMlevelRecursiveBisection
#define MCMlevelEdgeBisection		Metis_MCMlevelEdgeBisection
#define MCHMlevelEdgeBisection		Metis_MCHMlevelEdgeBisection


/* mrefine.c */
#define MocRefine2Way			Metis_MocRefine2Way
#define MocAllocate2WayPartitionMemory	Metis_MocAllocate2WayPartitionMemory
#define MocCompute2WayPartitionParams	Metis_MocCompute2WayPartitionParams
#define MocProject2WayPartition		Metis_MocProject2WayPartition


/* mrefine2.c */
#define MocRefine2Way2			Metis_MocRefine2Way2


/* mutil.c */
#define AreAllVwgtsBelow		Metis_AreAllVwgtsBelow
#define AreAnyVwgtsBelow		Metis_AreAnyVwgtsBelow
#define AreAllVwgtsAbove		Metis_AreAllVwgtsAbove
#define ComputeLoadImbalance		Metis_ComputeLoadImbalance
#define AreAllBelow			Metis_AreAllBelow


/* myqsort.c */
#define iidxsort			Metis_iidxsort
#define iintsort			Metis_iintsort
#define ikeysort			Metis_ikeysort
#define ikeyvalsort			Metis_ikeyvalsort


/* ometis.c */
#define MlevelNestedDissection		Metis_MlevelNestedDissection
#define MlevelNestedDissectionCC	Metis_MlevelNestedDissectionCC
#define MlevelNodeBisectionMultiple	Metis_MlevelNodeBisectionMultiple
#define MlevelNodeBisection		Metis_MlevelNodeBisection
#define SplitGraphOrder			Metis_SplitGraphOrder
#define MMDOrder			Metis_MMDOrder
#define SplitGraphOrderCC		Metis_SplitGraphOrderCC


/* parmetis.c */
#define MlevelNestedDissectionP		Metis_MlevelNestedDissectionP


/* pmetis.c */
#define MlevelRecursiveBisection	Metis_MlevelRecursiveBisection
#define MlevelEdgeBisection		Metis_MlevelEdgeBisection
#define SplitGraphPart			Metis_SplitGraphPart
#define SetUpSplitGraph			Metis_SetUpSplitGraph


/* pqueue.c */
#define PQueueInit			Metis_PQueueInit
#define PQueueReset			Metis_PQueueReset
#define PQueueFree			Metis_PQueueFree
#define PQueueInsert			Metis_PQueueInsert
#define PQueueDelete			Metis_PQueueDelete
#define PQueueUpdate			Metis_PQueueUpdate
#define PQueueUpdateUp			Metis_PQueueUpdateUp
#define PQueueGetMax			Metis_PQueueGetMax
#define PQueueSeeMax			Metis_PQueueSeeMax
#define CheckHeap			Metis_CheckHeap


/* refine.c */
#define Refine2Way			Metis_Refine2Way
#define Allocate2WayPartitionMemory	Metis_Allocate2WayPartitionMemory
#define Compute2WayPartitionParams	Metis_Compute2WayPartitionParams
#define Project2WayPartition		Metis_Project2WayPartition


/* separator.c */
#define ConstructSeparator		Metis_ConstructSeparator
#define ConstructMinCoverSeparator0	Metis_ConstructMinCoverSeparator0
#define ConstructMinCoverSeparator	Metis_ConstructMinCoverSeparator


/* sfm.c */
#define FM_2WayNodeRefine		Metis_FM_2WayNodeRefine
#define FM_2WayNodeRefineEqWgt		Metis_FM_2WayNodeRefineEqWgt
#define FM_2WayNodeRefine_OneSided	Metis_FM_2WayNodeRefine_OneSided
#define FM_2WayNodeBalance		Metis_FM_2WayNodeBalance
#define ComputeMaxNodeGain		Metis_ComputeMaxNodeGain


/* srefine.c */
#define Refine2WayNode			Metis_Refine2WayNode
#define Allocate2WayNodePartitionMemory	Metis_Allocate2WayNodePartitionMemory
#define Compute2WayNodePartitionParams	Metis_Compute2WayNodePartitionParams
#define Project2WayNodePartition	Metis_Project2WayNodePartition


/* stat.c */
#define ComputePartitionInfo		Metis_ComputePartitionInfo
#define ComputePartitionBalance		Metis_ComputePartitionBalance
#define ComputeElementBalance		Metis_ComputeElementBalance


/* subdomains.c */
#define Random_KWayEdgeRefineMConn	Metis_Random_KWayEdgeRefineMConn
#define Greedy_KWayEdgeBalanceMConn	Metis_Greedy_KWayEdgeBalanceMConn
#define PrintSubDomainGraph		Metis_PrintSubDomainGraph
#define ComputeSubDomainGraph		Metis_ComputeSubDomainGraph
#define EliminateSubDomainEdges		Metis_EliminateSubDomainEdges
#define MoveGroupMConn			Metis_MoveGroupMConn
#define EliminateComponents		Metis_EliminateComponents
#define MoveGroup			Metis_MoveGroup


/* timing.c */
#define InitTimers			Metis_InitTimers
#define PrintTimers			Metis_PrintTimers
#define seconds				Metis_seconds


/* util.c */
#define errexit				Metis_errexit
#define GKfree				Metis_GKfree
#ifndef DMALLOC
#define imalloc				Metis_imalloc
#define idxmalloc			Metis_idxmalloc
#define fmalloc				Metis_fmalloc
#define ismalloc			Metis_ismalloc
#define idxsmalloc			Metis_idxsmalloc
#define GKmalloc			Metis_GKmalloc
#endif
#define iset				Metis_iset
#define idxset				Metis_idxset
#define sset				Metis_sset
#define iamax				Metis_iamax
#define idxamax				Metis_idxamax
#define idxamax_strd			Metis_idxamax_strd
#define samax				Metis_samax
#define samax2				Metis_samax2
#define idxamin				Metis_idxamin
#define samin				Metis_samin
#define idxsum				Metis_idxsum
#define idxsum_strd			Metis_idxsum_strd
#define idxadd				Metis_idxadd
#define charsum				Metis_charsum
#define isum				Metis_isum
#define ssum				Metis_ssum
#define ssum_strd			Metis_ssum_strd
#define sscale				Metis_sscale
#define snorm2				Metis_snorm2
#define sdot				Metis_sdot
#define saxpy				Metis_saxpy
#define RandomPermute			Metis_RandomPermute
#define ispow2				Metis_ispow2
#define InitRandom			Metis_InitRandom
#define log2				Metis_log2





