
#module load e4s
#spack env activate gcc
#spack load cdo

hist_path='/pscratch/sd/h/huiwan/vis/ncic/cases/NUMLIQ_budget_0_F2010/tests/XS_1x1_nmonths/run/'
hist_path='/pscratch/sd/h/huiwan/vis/ncic/cases/NUMLIQ_budget_0_F2010/tests/XS_1x2_ndays/run/'
hist_file='NUMLIQ_budget_0_F2010.eam.h0.2009-10-02-00000.nc'
backhere=`pwd`

cd $hist_path

# Elements for constructing output variable names

cnd_prefix='cnd01'
qoi_name='NUMLIQ'
inc_suffix='_inc'
int_label=''
int_label='_v'

# Loop over all checkpoints get info using cdo

for chkpt in \
             'MCTCPL'   'PACINI'   'CHEM'      'AERDRYRM' \
             'GWDRAG'   'NDG'      'DYNEND'    'DEEPCU' \
             'CLDMAC01' 'CLDAER01' 'ACTDIAG01' 'CLDMIC01' \
             'CLDMAC02' 'CLDAER02' 'ACTDIAG02' 'CLDMIC02' \
             'CLDMAC03' 'CLDAER03' 'ACTDIAG03' 'CLDMIC03' \
             'CLDMAC04' 'CLDAER04' 'ACTDIAG04' 'CLDMIC04' \
             'CLDMAC05' 'CLDAER05' 'ACTDIAG05' 'CLDMIC05' \
             'CLDMAC06' 'CLDAER06' 'ACTDIAG06' 'CLDMIC06' \
             'AERWETRM' 'PBCDIAG' 
do

  echo -e '\n\n\n'
  varname=${cnd_prefix}'_'${qoi_name}${int_label}'_'${chkpt}${inc_suffix}
  cdo infov -selname,${varname} $hist_file
 #cdo infov -selname,${varname} -sellevel,0.868939 $hist_file
 #cdo -selname,${varname} -sellevel,0.868939 $hist_file $backhere/out.nc

done

cd $backhere
