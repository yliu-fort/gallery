/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2015 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "BoomerAMG.H"
#include "PstreamGlobals.H"
#include "processorLduInterface.H"
#include "processorGAMGInterfaceField.H"

// MPI library
#include <mpi.h>

// HYPRE libs
#include <math.h>

// static class
#include "Hypre_helper.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(BoomerAMG, 0);

    lduMatrix::solver::addsymMatrixConstructorToTable<BoomerAMG>
        addBoomerAMGSymMatrixConstructorToTable_;

}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::BoomerAMG::BoomerAMG
(
    const word& fieldName,
    const lduMatrix& matrix,
    const FieldField<Field, scalar>& interfaceBouCoeffs,
    const FieldField<Field, scalar>& interfaceIntCoeffs,
    const lduInterfaceFieldPtrsList& interfaces,
    const dictionary& solverControls
)
:
    lduMatrix::solver
    (
        fieldName,
        matrix,
        interfaceBouCoeffs,
        interfaceIntCoeffs,
        interfaces,
        solverControls
    )
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Hypre_helper hy;

Foam::solverPerformance Foam::BoomerAMG::solve
(
    scalarField& psi,
    const scalarField& source,
    const direction cmpt
) const
{

    //Info<< "   Greeting from BoomerAMG Solver..." << endl;
    //label numprocs = Pstream::nProcs(UPstream::worldComm);

    int myRank(0), numprocs(1);
    MPI_Comm mycomm(NULL);
    MPI_Status ierr;
    // Check if run in parallel...
    bool isParallel(PstreamGlobals::MPICommunicators_.size());
    //Info<< "   Check if run in parallel mode..." << (isParallel?"Yes":"No") << endl;

    if(isParallel)
    {
	// Retrieve MPI handle from lable
	mycomm = PstreamGlobals::MPICommunicators_[matrix().mesh().comm()];
	// Test the handle via requesting info

	MPI_Comm_size(mycomm, &numprocs);
	MPI_Comm_rank(mycomm, &myRank);

	//Pout<< "   n procs started:" << numprocs << endl;
    }else{}
    

    // Use printf instead of Info to output msg from all processors
    //Info<< "   MPIrun: I am " << myRank << " of " << numprocs << endl;
    //printf("I am %d of %d\n",myRank, numprocs);

    // --- Setup class containing solver performance data
    solverPerformance solverPerf
    (
        lduMatrix::preconditioner::getName(controlDict_) + typeName,
        fieldName_
    );

    // Ax = b, psi -> x, source -> b
    // res = b - Ax, rA -> res
    // Matrix_ -> A, need to retrieve spmatrix info from that class..
    label nCells = psi.size();

	if(hy.first)
	{
	// Compute global index offset...
	Info<< "   Calculate global index offset...";
	label global_offset_lower, global_offset_upper, my_tag_a(28);
	{
	    if(myRank == 0){
		global_offset_lower = 0;
		global_offset_upper = global_offset_lower + nCells;
		if(isParallel)
		    MPI_Ssend(&global_offset_upper, 1, MPI_INT, myRank+1, my_tag_a, mycomm); // to rank 1
	    }else{
		MPI_Recv(&global_offset_lower, 1, MPI_INT, myRank-1, my_tag_a, mycomm, &ierr); // Recv from lower rank
		// Compute offset...
		global_offset_upper = global_offset_lower + nCells;
		if(myRank < numprocs - 1) // Not the last one processor...
		    MPI_Ssend(&global_offset_upper, 1, MPI_INT, myRank+1, my_tag_a, mycomm);// Send computed offset to higher rank
	    }
	}
	Info << "ok" << endl;
	//Pout << "Global index offset: [" << global_offset_lower << ", " << global_offset_upper << "]" << endl;

	// Up to here we have retrieved owner and neighbour proc for all processor interfaces
	// next step will be convert the info to global cell index and merge the result into
	// derived ldu matrix.
	// To be able to do this, each processor should either send or receive info from another
	// to compute the global cell index. This requires to call some mpi related subroutines
	// GAMG should have done something similar...
	// see processorGAMGInterfaceField.C

	// Count the number of processor faces...
	Info<< "   Count the number of processor faces...";
	label count(0);
	    forAll(interfaces_, patchi)
	    {
		if (interfaces_.set(patchi))
		{
		    const labelUList& pa = matrix_.lduAddr().patchAddr(patchi);
		    count += pa.size();
		}
	}
	MPI_Barrier(mycomm);
	Info << "ok" << endl;
	// Retrive global index for processor faces on each processor...
	// Better new a container to store owner and values...
	//Pout << "Num of processor faces..." << count << endl;
	Info<< "   Exchange processor interface infomation...";
	label _neighbProcGlobalFaceIndexBuf[count], 
	_myProcGlobalFaceIndexBuf[count], 
	my_tag_b(100);
	label* _neighbProcGlobalFaceIndexBufPtr = &_neighbProcGlobalFaceIndexBuf[0];
	label* _myProcGlobalFaceIndexBufPtr = &_myProcGlobalFaceIndexBuf[0];
	scalar _ProcGlobalFaceCoeffBuf[count];
	scalar* _ProcGlobalFaceCoeffBufPtr = &_ProcGlobalFaceCoeffBuf[0];

	    forAll(interfaces_, patchi)
	    {
		if (interfaces_.set(patchi))
		{
		    const processorLduInterfaceField& p =
			 refCast<const processorLduInterfaceField>(interfaces_[patchi]);

		    const label neighbRank = p.neighbProcNo();

		    //word thisInterfaceType = interfaces_[patchi].interfaceFieldType();
		    const labelUList& pa = matrix_.lduAddr().patchAddr(patchi);
		    const scalarField& pCoeffs = interfaceBouCoeffs_[patchi];

		    // How many faces in this patch...
		    int count_ = pa.size();
		    
		    // Init a new list and copy pa element-wise-> add global offset
		    // Notice: this includes an assumption that everything(include sublists) is sorted.

		    forAll(pa, elem){
			*(_myProcGlobalFaceIndexBufPtr+elem) = pa[elem] + global_offset_lower;
			*(_ProcGlobalFaceCoeffBufPtr+elem) = -pCoeffs[elem];
		    }

		////Pout << " Global pa... " << count_ <<", "<<global_pa << endl;

		    // Send to higher rank first, then receive from lower rank...
		    // This is a lazy condition... might need to modify for cyclic boundaries.
		    if(myRank > neighbRank){
			MPI_Recv(_neighbProcGlobalFaceIndexBufPtr, count_, MPI_INT, neighbRank, my_tag_b+neighbRank, mycomm, &ierr);
			MPI_Ssend(_myProcGlobalFaceIndexBufPtr, count_, MPI_INT, neighbRank, my_tag_b+myRank, mycomm);
		    }else{
			MPI_Ssend(_myProcGlobalFaceIndexBufPtr, count_, MPI_INT, neighbRank, my_tag_b+myRank, mycomm);
			MPI_Recv(_neighbProcGlobalFaceIndexBufPtr, count_, MPI_INT, neighbRank, my_tag_b+neighbRank, mycomm, &ierr);
		    }

		    // Increment the ptr
		    _neighbProcGlobalFaceIndexBufPtr += count_;
		    _myProcGlobalFaceIndexBufPtr += count_;
		    _ProcGlobalFaceCoeffBufPtr += count_;

		}
	    }

	Info << "ok" << endl;
	// Create a label list to handle data...
	labelUList neighbProcGlobalFaceIndex(_neighbProcGlobalFaceIndexBuf, count);
	labelUList myProcGlobalFaceIndex(_myProcGlobalFaceIndexBuf, count);
	scalarUList procInterfaceCoeff(_ProcGlobalFaceCoeffBuf, count);
	//Pout << "Global cell for processor faces..." << myProcGlobalFaceIndex << ", " <<neighbProcGlobalFaceIndex << endl;

	// Build HYPRE-PARCSR handle...
	Info << "   Construct Hypre IJMatrix handle...";
	hy.ilower = global_offset_lower;
	hy.iupper = global_offset_upper-1;
	hy.jlower = global_offset_lower;
	hy.jupper = global_offset_upper-1;

	hy.jlower = hy.jlower<neighbProcGlobalFaceIndex.first()?hy.jlower:neighbProcGlobalFaceIndex.first();
	hy.jupper = hy.jupper>neighbProcGlobalFaceIndex.last()?hy.jupper:neighbProcGlobalFaceIndex.last();

	//Pout << "Constructing PARCSR, hy.ilower=" << hy.ilower << " hy.iupper=" << hy.iupper << " hy.jlower=" << hy.jlower << " hy.jupper=" << hy.jupper << endl;

	   /* Create the matrix.
	      Note that this is a square matrix, so we indicate the row partition
	      size twice (since number of rows = number of cols) */
	   HYPRE_IJMatrixCreate(mycomm, hy.ilower, hy.iupper, hy.ilower, hy.iupper, &hy.A);

	   /* Choose a parallel csr format storage (see the User's Manual) */
	   HYPRE_IJMatrixSetObjectType(hy.A, HYPRE_PARCSR);

	   /* Initialize before setting coefficients */
	   HYPRE_IJMatrixInitialize(hy.A);

	{
	    const scalar* __restrict__ diagPtr = matrix_.diag().begin();
	    const scalar* __restrict__ lowerPtr = matrix_.lower().begin();
	    const scalar* __restrict__ upperPtr = matrix_.upper().begin();
	    const label* __restrict__ lPtr = matrix_.lduAddr().lowerAddr().begin();
	    const label* __restrict__ uPtr = matrix_.lduAddr().upperAddr().begin();

	    // Ptr for processor interfaces..
	    const label* procIPtr = myProcGlobalFaceIndex.begin();
	    const label* procJPtr = neighbProcGlobalFaceIndex.begin();
	    const scalar* procValPtr = procInterfaceCoeff.begin();

	    int count = 0;
		// Lower
		forAll(matrix_.lduAddr().upperAddr(), i)
		{
		    int row(uPtr[i] + global_offset_lower), col(lPtr[i] + global_offset_lower), ncols(1);
		    double val(lowerPtr[i]);

		    HYPRE_IJMatrixSetValues(hy.A, 1, &ncols, &row, &col, &val); // SpMatrix, nrows, ncols, rowPtr, colPtr, val
		    count++;
		}

		// Upper
		forAll(matrix_.lduAddr().lowerAddr(), i)
		{
		    int row(lPtr[i] + global_offset_lower), col(uPtr[i] + global_offset_lower), ncols(1);
		    double val(upperPtr[i]);

		    HYPRE_IJMatrixSetValues(hy.A, 1, &ncols, &row, &col, &val); // SpMatrix, nrows, ncols, rowPtr, colPtr, val
		    count++;
		}

		// Diag
		forAll(matrix_.diag(), i)
		{
		    int row(i + global_offset_lower), col(i + global_offset_lower), ncols(1);
		    double val(diagPtr[i]);

		    HYPRE_IJMatrixSetValues(hy.A, 1, &ncols, &row, &col, &val); // SpMatrix, nrows, ncols, rowPtr, colPtr, val
		    count++;
		}

		// Processor Interfaces
		forAll(myProcGlobalFaceIndex, i)
		{
		    int row(procIPtr[i]), col(procJPtr[i]), ncols(1);
		    double val(procValPtr[i]);

		    HYPRE_IJMatrixSetValues(hy.A, 1, &ncols, &row, &col, &val); // SpMatrix, nrows, ncols, rowPtr, colPtr, val
		    count++;
		}

	}
	//MPI_Barrier(mycomm);
	Info << "ok" << endl;
	// Bind to HYPRE interface...
	Info << "   Bind handle to Hypre solver interface...";

	   /* Assemble after setting the coefficients */
	   HYPRE_IJMatrixAssemble(hy.A);

	   /* Get the parcsr matrix object to use */
	   HYPRE_IJMatrixGetObject(hy.A, (void**) &hy.parcsr_A);

	   // Bind solver handle...
	   HYPRE_BoomerAMGCreate(&hy.solver);

	}

   /* Create the rhs and solution */
   HYPRE_IJVector b;
   HYPRE_ParVector par_b;
   HYPRE_IJVector x;
   HYPRE_ParVector par_x;
   HYPRE_IJVectorCreate(mycomm, hy.ilower, hy.iupper,&b);
   HYPRE_IJVectorSetObjectType(b, HYPRE_PARCSR);
   HYPRE_IJVectorInitialize(b);

   HYPRE_IJVectorCreate(mycomm, hy.ilower, hy.iupper,&x);
   HYPRE_IJVectorSetObjectType(x, HYPRE_PARCSR);
   HYPRE_IJVectorInitialize(x);

   /* Set the rhs values and the solution */

{
   int rows[nCells];
   for(int i = 0;i < nCells; i++)
    {
	rows[i] = i + hy.ilower;
    }

   HYPRE_IJVectorSetValues(b, nCells, rows, &source[0]);
   HYPRE_IJVectorSetValues(x, nCells, rows, &psi[0]);
}

   HYPRE_IJVectorAssemble(b);
   HYPRE_IJVectorGetObject(b, (void **) &par_b);

   HYPRE_IJVectorAssemble(x);
   HYPRE_IJVectorGetObject(x, (void **) &par_x);

// Solve matrix...
{
      int num_iterations;
      double final_res_norm;

      /* Create solver */
      //HYPRE_BoomerAMGCreate(&solver);
	if(hy.first)
	{
	      /* Set some parameters (See Reference Manual for more parameters) */
	      HYPRE_BoomerAMGSetPrintLevel(hy.solver, 1);  /* print solve info + parameters */
	      HYPRE_BoomerAMGSetOldDefault(hy.solver); /* Falgout coarsening with modified classical interpolaiton */
	      HYPRE_BoomerAMGSetRelaxType(hy.solver, 3);   /* G-S/Jacobi hybrid relaxation */
	      HYPRE_BoomerAMGSetRelaxOrder(hy.solver, 1);   /* uses C/F relaxation */
	      HYPRE_BoomerAMGSetNumSweeps(hy.solver, 1);   /* Sweeeps on each level */
	      HYPRE_BoomerAMGSetMaxLevels(hy.solver, 20);  /* maximum number of levels */
	      //HYPRE_BoomerAMGSetStrongThreshold(hy.solver, 0.5); /* affect complexity of the solver*/
	      HYPRE_BoomerAMGSetTol(hy.solver, tolerance_);      /* conv. tolerance */

	      /* Now setup and solve! */
	      HYPRE_BoomerAMGSetup(hy.solver, hy.parcsr_A, par_b, par_x);

	     // set first = false
	     hy.first = false;
	}
      HYPRE_BoomerAMGSolve(hy.solver, hy.parcsr_A, par_b, par_x);

      /* Run info - needed logging turned on */
      HYPRE_BoomerAMGGetNumIterations(hy.solver, &num_iterations);
      HYPRE_BoomerAMGGetFinalRelativeResidualNorm(hy.solver, &final_res_norm);

      /* Destroy solver */
      //HYPRE_BoomerAMGDestroy(solver);

   // Retrieve solution...
   solverPerf.finalResidual() = final_res_norm;
   solverPerf.nIterations() = num_iterations;

   int rows[nCells];
   
   for(int i = 0;i < nCells; i++)
    {
	rows[i] = i + hy.ilower;
    }

   HYPRE_IJVectorGetValues(x, nCells, rows, &psi[0]);

}

//scalarUList psi_buf(_psi, nCells);
////Pout << "Psi " << psi_buf << endl;


   /* Clean up */
   //HYPRE_IJMatrixDestroy(A);
   HYPRE_IJVectorDestroy(b);
   HYPRE_IJVectorDestroy(x);

    return solverPerf;
}


// ************************************************************************* //
