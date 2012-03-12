#include"ccsd.h"
#include"blas.h"
#ifdef _OPENMP
   #include<omp.h>
#endif

using namespace psi;

namespace psi{
PsiReturnType triples(boost::shared_ptr<psi::CoupledCluster>ccsd,Options&options);

PsiReturnType triples(boost::shared_ptr<psi::CoupledCluster>ccsd,Options&options){

  fprintf(outfile,"\n");
  fprintf(outfile, "        *******************************************************\n");
  fprintf(outfile, "        *                                                     *\n");
  fprintf(outfile, "        *                      CCSD(T)                        *\n");
  fprintf(outfile, "        *                                                     *\n");
  fprintf(outfile, "        *******************************************************\n");
  fprintf(outfile,"\n");
  fflush(outfile);

  int o = ccsd->ndoccact;
  int v = ccsd->nvirt_no;

  double **Rii = ccsd->wfn_->Rii->pointer();
  double *t1 = ccsd->t1;
  double *F  = ccsd->eps;
  double *E2ijak,**E2abci;
  E2ijak = (double*)malloc(o*o*o*v*sizeof(double));
  int nthreads = 1;
  #ifdef _OPENMP
      nthreads = omp_get_max_threads();
  #endif

  if (options["NUM_THREADS"].has_changed())
     nthreads = options.get_int("NUM_THREADS");

  long int memory = Process::environment.get_memory();
  if (options["CCMEMORY"].has_changed()){
     memory  = options.get_int("CCMEMORY");
     memory *= (long int)1024*1024;
  }
  memory -= 8L*(2L*o*o*v*v+o*o*o*v+o*v+3L*nthreads*v*v*v);

  fprintf(outfile,"        num_threads =             %9i\n",nthreads);
  fprintf(outfile,"        available memory =     %9.2lf mb\n",memory/1024./1024.);
  fprintf(outfile,"        memory requirements =  %9.2lf mb\n",
           8.*(2.*o*o*v*v+1.*o*o*o*v+(3.*nthreads)*v*v*v+1.*o*v)/1024./1024.);
  fprintf(outfile,"\n");
  fflush(outfile);

  bool threaded = true;
  if (memory<0){
     memory += (nthreads-1)*8L*3L*v*v*v;
     if (nthreads==1){
        fprintf(outfile,"        Error: not enough memory.\n");
        fprintf(outfile,"\n");
        fprintf(outfile,"        (T) requires at least %7.2lf mb\n",
             8.*(2.*o*o*v*v+1.*o*o*o*v+3.*v*v*v+1.*o*v)/1024./1024.);
        fprintf(outfile,"\n");
        fflush(outfile);
        return Failure;
     }
     threaded = false;
     nthreads = 1;
     fprintf(outfile,"        Not enough memory for explicit threading ... \n");
     fprintf(outfile,"\n");
     fprintf(outfile,"        memory requirements =  %9.2lf mb\n",
              8.*(2.*o*o*v*v+1.*o*o*o*v+(3.)*v*v*v+1.*o*v)/1024./1024.);
     fprintf(outfile,"\n");
     fflush(outfile);
  }

  E2abci = (double**)malloc(nthreads*sizeof(double*));
  // some v^3 intermediates
  double **Z  = (double**)malloc(nthreads*sizeof(double*));
  double **Z2 = (double**)malloc(nthreads*sizeof(double*));

  boost::shared_ptr<PSIO> psio(new PSIO());
  double*tempE2=(double*)malloc(o*o*o*v*sizeof(double));
  psio->open(PSIF_IJAK,PSIO_OPEN_OLD);
  psio->read_entry(PSIF_IJAK,"E2ijak",(char*)&tempE2[0],o*o*o*v*sizeof(double));
  psio->close(PSIF_IJAK,1);
  for (int i=0; i<o*o*o; i++){
      for (int a=0; a<v; a++){
          E2ijak[a*o*o*o+i] = tempE2[i*v+a];
      }
  }
  free(tempE2);

  for (int i=0; i<nthreads; i++){
      E2abci[i] = (double*)malloc(v*v*v*sizeof(double));
      Z[i]      = (double*)malloc(v*v*v*sizeof(double));
      Z2[i]     = (double*)malloc(v*v*v*sizeof(double));
  }


  double *tempt = (double*)malloc(o*o*v*v*sizeof(double));

  if (ccsd->t2_on_disk){
     ccsd->tb = (double*)malloc(o*o*v*v*sizeof(double));
     psio->open(PSIF_T2,PSIO_OPEN_OLD);
     psio->read_entry(PSIF_T2,"t2",(char*)&ccsd->tb[0],o*o*v*v*sizeof(double));
     psio->close(PSIF_T2,1);
  }

  //for (int a=0; a<v*v; a++){
  //    F_DCOPY(o*o,ccsd->tb+a*o*o,1,tempt+a,v*v);
  //}
  F_DCOPY(o*o*v*v,ccsd->tb,1,tempt,1);

  // might as well use t2's memory
  double*E2klcd = ccsd->tb;
  psio->open(PSIF_KLCD,PSIO_OPEN_OLD);
  psio->read_entry(PSIF_KLCD,"E2klcd", (char*)&E2klcd[0],o*o*v*v*sizeof(double));
  psio->close(PSIF_KLCD,1);

  double *etrip = (double*)malloc(nthreads*sizeof(double));
  for (int i=0; i<nthreads; i++) etrip[i] = 0.0;
  fprintf(outfile,"        Computing (T) correction...\n");
  fprintf(outfile,"\n");
  fprintf(outfile,"        %% complete  total time\n");
  fflush(outfile);

  time_t stop,start = time(NULL);
  int pct10,pct20,pct30,pct40,pct50,pct60,pct70,pct80,pct90;
  pct10=pct20=pct30=pct40=pct50=pct60=pct70=pct80=pct90=0;

  int nabc = 0;
  for (int a=0; a<v; a++){
      for (int b=0; b<v; b++){
          for (int c=0; c<v; c++){
              nabc++;
          }
      }
  }
  int**abc = (int**)malloc(nabc*sizeof(int*));
  nabc = 0;
  for (int a=0; a<v; a++){
      for (int b=0; b<v; b++){
          for (int c=0; c<v; c++){
              abc[nabc] = (int*)malloc(3*sizeof(int));
              abc[nabc][0] = a;
              abc[nabc][1] = b;
              abc[nabc][2] = c;
              nabc++;
          }
      }
  }
  fprintf(outfile,"        Number of abc combinations: %i\n",nabc);
  fprintf(outfile,"\n");
  fflush(outfile);
  for (int i=0; i<nthreads; i++) etrip[i] = 0.0;

  /**
    *  if there is enough memory to explicitly thread, do so
    */
  if (threaded){
     #pragma omp parallel for schedule (dynamic) num_threads(nthreads)
     for (int ind=0; ind<nabc; ind++){
         int a = abc[ind][0];
         int b = abc[ind][1];
         int c = abc[ind][2];

         int thread = 0;
         #ifdef _OPENMP
             thread = omp_get_thread_num();
         #endif

         boost::shared_ptr<PSIO> mypsio(new PSIO());
         mypsio->open(PSIF_ABCI4,PSIO_OPEN_OLD);
         psio_address addr = psio_get_address(PSIO_ZERO,(long int)(b*v*v*o+c*v*o)*sizeof(double));
         mypsio->read(PSIF_ABCI4,"E2abci4",(char*)&E2abci[thread][0],o*v*sizeof(double),addr,&addr);
        
         // (1)
         F_DGEMM('t','t',o,o*o,v,1.0,E2abci[thread],v,tempt+a*o*o*v,o*o,0.0,Z[thread],o); 
         // (ikj)(acb)
         F_DGEMM('t','n',o,o*o,o,-1.0,tempt+c*o*o*v+a*o*o,o,E2ijak+b*o*o*o,o,1.0,Z[thread],o); 
         //for (int i=0; i<o; i++){
         //    for (int j=0; j<o; j++){
         //        for (int k=0; k<o; k++){
         //            double dum = 0.0;
                     //for (int e=0; e<v; e++){
                     //    dum += E2abci[thread][k*v+e] * tempt[a*o*o*v+e*o*o+i*o+j];
                     //}
                     //for (int m=0; m<o; m++){
                     //    dum -= E2ijak[b*o*o*o+i*o*o+j*o+m] * tempt[c*o*o*v+a*o*o+k*o+m];
                     //}
                     //Z[thread][i*o*o+j*o+k] += dum;
         //        }
         //    }
         //}

         addr = psio_get_address(PSIO_ZERO,(long int)(a*v*v*o+c*v*o)*sizeof(double));
         mypsio->read(PSIF_ABCI4,"E2abci4",(char*)&E2abci[thread][0],o*v*sizeof(double),addr,&addr);
         //(ab)(ij)
         F_DGEMM('t','t',o,o*o,v,1.0,E2abci[thread],v,tempt+b*o*o*v,o*o,0.0,Z2[thread],o);
         //(ab)(ij)
         F_DGEMM('t','n',o*o,o,o,-1.0,E2ijak+c*o*o*o,o,tempt+b*o*o*v+a*o*o,o,1.0,Z2[thread],o*o);
         for (int i=0; i<o; i++){
             for (int j=0; j<o; j++){
                 F_DAXPY(o,1.0,Z2[thread]+j*o*o+i*o,1,Z[thread]+i*o*o+j*o,1);
                 //for (int k=0; k<o; k++){
                     //double dum = 0.0;
                     //for (int e=0; e<v; e++){
                     //    dum += E2abci[thread][k*v+e] * tempt[b*o*o*v+e*o*o+j*o+i];
                     //}
                     //for (int m=0; m<o; m++){
                     //    dum -= E2ijak[c*o*o*o+i*o*o+k*o+m] * tempt[b*o*o*v+a*o*o+j*o+m];
                     //}
                     //Z[thread][i*o*o+j*o+k] += Z2[thread][j*o*o+i*o+k];
                 //}
             }
         }

         addr = psio_get_address(PSIO_ZERO,(long int)(c*v*v*o+b*v*o)*sizeof(double));
         mypsio->read(PSIF_ABCI4,"E2abci4",(char*)&E2abci[thread][0],o*v*sizeof(double),addr,&addr);
         //(bc)(jk)
         F_DGEMM('t','t',o,o*o,v,1.0,E2abci[thread],v,tempt+a*o*o*v,o*o,0.0,Z2[thread],o);
         //(bc)(jk)
         F_DGEMM('t','n',o*o,o,o,-1.0,E2ijak+b*o*o*o,o,tempt+a*o*o*v+c*o*o,o,1.0,Z2[thread],o*o);
         for (int i=0; i<o; i++){
             for (int j=0; j<o; j++){
                 F_DAXPY(o,1.0,Z2[thread]+i*o*o+j,o,Z[thread]+i*o*o+j*o,1);
                // for (int k=0; k<o; k++){
                //     double dum = 0.0;
                //     for (int e=0; e<v; e++){
                //         dum += E2abci[thread][j*v+e] * tempt[a*o*o*v+e*o*o+i*o+k];
                //     }
                //     for (int m=0; m<o; m++){
                //         dum -= E2ijak[b*o*o*o+k*o*o+j*o+m] * tempt[a*o*o*v+c*o*o+i*o+m];
                //     }
                //    Z[thread][i*o*o+j*o+k] += dum;
                // }
             }
         }
         addr = psio_get_address(PSIO_ZERO,(long int)(b*v*v*o+a*v*o)*sizeof(double));
         mypsio->read(PSIF_ABCI4,"E2abci4",(char*)&E2abci[thread][0],o*v*sizeof(double),addr,&addr);
         //(ac)(ik)
         F_DGEMM('t','t',o,o*o,v,1.0,E2abci[thread],v,tempt+c*o*o*v,o*o,0.0,Z2[thread],o);
         //(ac)(ik)
         F_DGEMM('t','n',o*o,o,o,-1.0,E2ijak+a*o*o*o,o,tempt+c*o*o*v+b*o*o,o,1.0,Z2[thread],o*o);
         //(1)
         F_DGEMM('t','t',o,o*o,o,-1.0,tempt+a*o*o*v+b*o*o,o,E2ijak+c*o*o*o,o*o,1.0,Z2[thread],o);
         for (int i=0; i<o; i++){
             for (int j=0; j<o; j++){
                 for (int k=0; k<o; k++){
                     double dum = 0.0;
                     //for (int e=0; e<v; e++){
                     //    dum += E2abci[thread][i*v+e] * tempt[c*o*o*v+e*o*o+k*o+j];
                     //}
                     //for (int m=0; m<o; m++){
                     //    dum -= E2ijak[a*o*o*o+j*o*o+i*o+m] * tempt[c*o*o*v+b*o*o+k*o+m];
                     //}
                     //  (1)
                     //for (int m=0; m<o; m++){
                     //    dum -= E2ijak[c*o*o*o+m*o*o+k*o+j] * tempt[a*o*o*v+b*o*o+i*o+m];
                     //}
                     Z[thread][i*o*o+j*o+k] += Z2[thread][k*o*o+j*o+i];
                 }
             }
         }
         addr = psio_get_address(PSIO_ZERO,(long int)(c*v*v*o+a*v*o)*sizeof(double));
         mypsio->read(PSIF_ABCI4,"E2abci4",(char*)&E2abci[thread][0],o*v*sizeof(double),addr,&addr);
         //(ijk)(abc)
         F_DGEMM('t','t',o,o*o,v,1.0,E2abci[thread],v,tempt+b*o*o*v,o*o,0.0,Z2[thread],o);
         F_DGEMM('t','n',o*o,o,o,-1.0,E2ijak+a*o*o*o,o,tempt+b*o*o*v+c*o*o,o,1.0,Z2[thread],o*o);
         //(ijk)(abc)
         for (int i=0; i<o; i++){
             for (int j=0; j<o; j++){
                 for (int k=0; k<o; k++){
                     double dum = 0.0;
                     //for (int e=0; e<v; e++){
                     //    dum += E2abci[thread][i*v+e] * tempt[b*o*o*v+e*o*o+j*o+k];
                     //}
                     //for (int m=0; m<o; m++){
                     //    dum -= E2ijak[a*o*o*o+k*o*o+i*o+m] * tempt[b*o*o*v+c*o*o+j*o+m];
                     //}
                     //Z[thread][i*o*o+j*o+k] += Z2[thread][j*o*o+k*o+i];
                 }
             }
         }
         //(ikj)(acb)
         addr = psio_get_address(PSIO_ZERO,(long int)(a*v*v*o+b*v*o)*sizeof(double));
         mypsio->read(PSIF_ABCI4,"E2abci4",(char*)&E2abci[thread][0],o*v*sizeof(double),addr,&addr);
         F_DGEMM('n','n',o*o,o,v,1.0,tempt+c*o*o*v,o*o,E2abci[thread],v,1.0,Z2[thread],o*o);
         for (int i=0; i<o; i++){
             for (int j=0; j<o; j++){
                 for (int k=0; k<o; k++){
                     //double dum = 0.0;
                     //for (int e=0; e<v; e++){
                     //    dum += E2abci[thread][j*v+e] * tempt[c*o*o*v+e*o*o+k*o+i];
                     //}
                     //Z[thread][i*o*o+j*o+k] += dum;
                     Z[thread][i*o*o+j*o+k] += Z2[thread][j*o*o+k*o+i];
                 }
             }
         }

//heyheyhey
         for (int i=0; i<o; i++){
             for (int j=0; j<o; j++){
                 for (int k=0; k<o; k++){
                     Z2[thread][i*o*o+j*o+k] = 4.0/3.0*Z[thread][i*o*o+j*o+k]
                                             -     2.0*Z[thread][i*o*o+k*o+j]
                                             + 2.0/3.0*Z[thread][j*o*o+k*o+i];
                 }
             }
         }
         for (int i=0; i<o; i++){
             for (int j=0; j<o; j++){
                 for (int k=0; k<o; k++){
                     double denom = F[i]+F[j]+F[k]-F[a+o]-F[b+o]-F[c+o];
                     Z[thread][i*o*o+j*o+k] /= denom;
                 }
             }
         }
         // transform Z and Z2 back to lmo basis
         F_DGEMM('n','t',o*o,o,o,1.0,Z[thread],o*o,&Rii[0][0],o,0.0,E2abci[thread],o*o);
         F_DCOPY(o*o*o,E2abci[thread],1,Z[thread],1);
         F_DGEMM('n','t',o*o,o,o,1.0,Z2[thread],o*o,&Rii[0][0],o,0.0,E2abci[thread],o*o);
         F_DCOPY(o*o*o,E2abci[thread],1,Z2[thread],1);
         for (int i=0; i<o; i++){
             for (int j=0; j<o; j++){
                 for (int k=0; k<o; k++){
                     etrip[thread] += Z[thread][i*o*o+j*o+k]*Z2[thread][i*o*o+j*o+k]*ccsd->wfn_->centralfac[i];
                 }
             }
         }
         
         // the (t) part ...
         for (int i=0; i<o; i++){
             double tai = t1[a*o+i];
             for (int j=0; j<o; j++){
                 double tbj = t1[b*o+j];
                 double E2iajb = E2klcd[i*v*v*o+a*v*o+j*v+b];
                 for (int k=0; k<o; k++){
                     Z2[thread][i*o*o+j*o+k] = (tai      *E2klcd[j*v*v*o+b*v*o+k*v+c] +
                                                tbj      *E2klcd[i*v*v*o+a*v*o+k*v+c] +
                                                t1[c*o+k]*E2iajb);
                 }
             }
         }
         for (int i=0; i<o; i++){
             for (int j=0; j<o; j++){
                 for (int k=0; k<o; k++){
                     E2abci[thread][i*o*o+j*o+k] = 4.0/3.0*Z2[thread][i*o*o+j*o+k]
                                                 -     2.0*Z2[thread][i*o*o+k*o+j]
                                                 + 2.0/3.0*Z2[thread][j*o*o+k*o+i];
                 }
             }
         }
         // transform Z2 back to lmo basis
         F_DGEMM('n','t',o*o,o,o,1.0,E2abci[thread],o*o,&Rii[0][0],o,0.0,Z2[thread],o*o);
         for (int i=0; i<o; i++){
             for (int j=0; j<o; j++){
                 for (int k=0; k<o; k++){
                     etrip[thread] += Z[thread][i*o*o+j*o+k]*Z2[thread][i*o*o+j*o+k]*ccsd->wfn_->centralfac[i];
                 }
             }
         }
         /*F_DCOPY(o*o*o,Z[thread],1,Z2[thread],1);
         for (int i=0; i<o; i++){
             double tai = t1[a*o+i];
             for (int j=0; j<o; j++){
                 int ij = 1+(i==j);
                 double tbj = t1[b*o+j];
                 double E2iajb = E2klcd[i*v*v*o+a*v*o+j*v+b];
                 for (int k=0; k<o; k++){
                     Z2[thread][i*o*o+j*o+k] += (tai      *E2klcd[j*v*v*o+b*v*o+k*v+c] +
                                                 tbj      *E2klcd[i*v*v*o+a*v*o+k*v+c] +
                                                 t1[c*o+k]*E2iajb);
                     Z2[thread][i*o*o+j*o+k] /= (ij + (j==k) + (i==k));
                 }
             }
         }

         for (int i=0; i<o; i++){
             for (int j=0; j<o; j++){
                 for (int k=0; k<o; k++){
                     long int ijk = i*o*o+j*o+k;
                     long int jik = j*o*o+i*o+k;
                     long int ikj = i*o*o+k*o+j;
                     long int kji = k*o*o+j*o+i;

                     E2abci[thread][ijk] = Z2[thread][ikj] + Z2[thread][jik] + Z2[thread][kji];
                 }
             }
         }
         double dabc = -F[a+o]-F[b+o]-F[c+o];
         // separate out these bits to save v^3 storage
         double tripval = 0.0;
         int abcfac = ( 2-((a==b)+(b==c)+(a==c)) );
         for (int i=0; i<o; i++){
             double dabci = dabc+F[i];
             for (int j=0; j<=i; j++){
                 double dabcij = dabci+F[j];
                 for (int k=0; k<=j; k++){

                     long int ijk = i*o*o+j*o+k;
                     long int jki = j*o*o+k*o+i;
                     long int kij = k*o*o+i*o+j;
                     long int ikj = i*o*o+k*o+j;
                     long int jik = j*o*o+i*o+k;
                     long int kji = k*o*o+j*o+i;
                     double dum      = Z[thread][ijk]*Z2[thread][ijk] + Z[thread][ikj]*Z2[thread][ikj]
                                     + Z[thread][jik]*Z2[thread][jik] + Z[thread][jki]*Z2[thread][jki]
                                     + Z[thread][kij]*Z2[thread][kij] + Z[thread][kji]*Z2[thread][kji];

                     dum            =  (E2abci[thread][ijk])
                                    * ((Z[thread][ijk] + Z[thread][jki] + Z[thread][kij])*-2.0
                                    +  (Z[thread][ikj] + Z[thread][jik] + Z[thread][kji]))
                                    + 3.0*dum;
                     double denom = dabcij+F[k];
                     tripval += dum/denom;
                 }
             }
         }
         etrip[thread] += tripval*abcfac;
         // the second bit
         for (int i=0; i<o; i++){
             for (int j=0; j<o; j++){
                 for (int k=0; k<o; k++){
                     long int ijk = i*o*o+j*o+k;
                     long int jki = j*o*o+k*o+i;
                     long int kij = k*o*o+i*o+j;

                     E2abci[thread][ijk]  = Z2[thread][ijk] + Z2[thread][jki] + Z2[thread][kij];
                 }
             }
         }
         tripval = 0.0;
         for (int i=0; i<o; i++){
             double dabci = dabc+F[i];
             for (int j=0; j<=i; j++){
                 double dabcij = dabci+F[j];
                 for (int k=0; k<=j; k++){
                     long int ijk = i*o*o+j*o+k;
                     long int jki = j*o*o+k*o+i;
                     long int kij = k*o*o+i*o+j;
                     long int ikj = i*o*o+k*o+j;
                     long int jik = j*o*o+i*o+k;
                     long int kji = k*o*o+j*o+i;

                     double dum     = (E2abci[thread][ijk])
                                    * (Z[thread][ijk] + Z[thread][jki] + Z[thread][kij]
                                    + (Z[thread][ikj] + Z[thread][jik] + Z[thread][kji])*-2.0);

                     double denom = dabcij+F[k];
                     tripval += dum/denom;
                 }
             }
         }
         etrip[thread] += tripval*abcfac;*/
         // print out update 
         if (thread==0){
            int print = 0;
            stop = time(NULL);
            if ((double)ind/nabc >= 0.1 && !pct10){      pct10 = 1; print=1;}
            else if ((double)ind/nabc >= 0.2 && !pct20){ pct20 = 1; print=1;}
            else if ((double)ind/nabc >= 0.3 && !pct30){ pct30 = 1; print=1;}
            else if ((double)ind/nabc >= 0.4 && !pct40){ pct40 = 1; print=1;}
            else if ((double)ind/nabc >= 0.5 && !pct50){ pct50 = 1; print=1;}
            else if ((double)ind/nabc >= 0.6 && !pct60){ pct60 = 1; print=1;}
            else if ((double)ind/nabc >= 0.7 && !pct70){ pct70 = 1; print=1;}
            else if ((double)ind/nabc >= 0.8 && !pct80){ pct80 = 1; print=1;}
            else if ((double)ind/nabc >= 0.9 && !pct90){ pct90 = 1; print=1;}
            if (print){
               fprintf(outfile,"              %3.1lf  %8d s\n",100.0*ind/nabc,(int)stop-(int)start);
               fflush(outfile);
            }
         }
         mypsio->close(PSIF_ABCI4,1);
         mypsio.reset();
     }
  }
  else{
     fprintf(outfile,"on the to do pile!\n");
     return Failure;
  }


  double et = 0.0;
  for (int i=0; i<nthreads; i++) et += etrip[i];

  fprintf(outfile,"\n");
  if (ccsd->scale_t == 1.0){
     fprintf(outfile,"        (T) energy                   %20.12lf\n",et);
  }
  else{
     fprintf(outfile,"                                                 unscaled               scaled\n");
     fprintf(outfile,"        (T) energy                   %20.12lf %20.12lf\n",et,et*ccsd->scale_t);
  }
  fprintf(outfile,"\n");
  if (ccsd->scale_t == 1.0)
     fprintf(outfile,"        CCSD(T) correlation energy   %20.12lf\n",ccsd->eccsd+et);
  else{
     fprintf(outfile,"                                                 unscaled               scaled\n");
     fprintf(outfile,"        CCSD(T) correlation energy   %20.12lf %20.12lf\n",ccsd->eccsd+et,ccsd->eccsd+et*ccsd->scale_t);
  }
  if (ccsd->scale_t == 1.0)
     fprintf(outfile,"      * CCSD(T) total energy         %20.12lf\n",ccsd->eccsd+et+ccsd->escf);
  else{
     fprintf(outfile,"                                                 unscaled               scaled\n");
     fprintf(outfile,"      * CCSD(T) total energy         %20.12lf %20.12lf\n",ccsd->eccsd+et+ccsd->escf,ccsd->eccsd+et*ccsd->scale_t+ccsd->escf);
  }
  fflush(outfile);
  ccsd->et = et;

  // free memory:
  if (ccsd->t2_on_disk){
     free(ccsd->tb);
  }
  free(E2ijak);
  for (int i=0; i<nthreads; i++){  
      free(E2abci[i]);
      free(Z[i]);
      free(Z2[i]);
  }
  free(Z);
  free(Z2);
  free(E2abci);
  free(etrip);
            
  return Success;
}


} // end of namespace



