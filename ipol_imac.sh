#!/bin/bash

echo

realpath() {
  echo "$(cd "$(dirname "$1")" && pwd)/$(basename "$1")"
}

ipoldir="$(realpath )"

starg=$1
enarg=$2
raddir="$(realpath $3)"
tempdir="$(realpath $4)"
doppdir=$5
declare -a simdirs=(${@:6})

ptype=mp4

mkdir -p $raddir $tempdir

declare -a mpopts=( "mp06" "mp08" "mp10" "mp16" "mp51" )
declare -a mpnames=( "wsm6" "thom" "morr" "wdm6" "p3" )

st=$(echo $starg | sed 's/[^0-9]*//g')
stdate=${st:0:8}
sttime=${st:8:4}
en=$(echo $enarg | sed 's/[^0-9]*//g')
endate=${en:0:8}
entime=${en:8:4}
stt=${stdate}_${sttime}
edt=${endate}_${entime}

configdir=$ipoldir/configtxt/${stt}_${edt}
mkdir -p $configdir

echo Determining data agency and type...
echo
sleep 3

station=$(basename $raddir)
if [[ "$station" == "CASAG" ]]; then
    agency='cwr'
elif [[ "$station" == "NPOL" ]]; then
    agency='olympex'
else
    agency='nexrad'
fi
echo $agency

if [[ "$(ls $raddir/* | head -n 1 | xargs basename)" == "wrfout"* ]]; then
	data='wrf'
    mp=$(ls $raddir/* | head -n 1 | xargs basename | cut -d '_' -f2)
    for ((ii=0;ii<${#mpopts[@]};ii++)); do
       [[ "${mpopts[ii]}" = "$mp" ]] && break
    done
    inputfile=input_${mpnames[ii]}_${stt}_${edt}.txt
    configfile=config_${mpnames[ii]}_${stt}_${edt}.txt
    echo $data
else
    if [ -z $simdirs ]; then
        data='obs'
        inputfile=input_${data}_${stt}_${edt}.txt
        configfile=config_${data}_${stt}_${edt}.txt
        echo $data
    else
        inputfile=input_obs_${stt}_${edt}.txt
        configfile=config_obs_${stt}_${edt}.txt
        echo "obs vs wrf"
    fi
fi

path2specs=$(realpath $raddir/../../../radar_specs/${station}_specs.txt)
#path2specs=$(realpath $raddir/../../${station}_specs.txt)

latcen=$(cat $path2specs | grep "latitude =" | cut -d '=' -f2 | cut -d ';' -f1 | xargs)
loncen=$(cat $path2specs | grep "longitude =" | cut -d '=' -f2 | cut -d ';' -f1 | xargs)
maxelev=$(cat $path2specs | grep "max_elev_angle =" | cut -d '=' -f2 | cut -d ';' -f1 | xargs)

echo
echo Selecting radar files for analysis in range $stt to $edt...
echo
sleep 3

tfile=input_${stt}_${edt}.txt

for filepath in $(ls $raddir/* | sort); do
    file=$(basename $filepath)
    if [[ "$(ls $raddir/* | head -n 1 | xargs basename)" == "wrfout"* ]]; then
        filedt=$(echo $file | cut -d '_' -f4 | tr -d '-')$(echo $file | cut -d '_' -f5 | cut -d '.' -f1 | tr -d ':')
    else
        filedt=$(echo $file | cut -d '_' -f2)$(echo $file | cut -d '_' -f3)
    fi
    if [ "$filedt" -ge "$(echo $stt | tr -d '_')00" ] && [ "$filedt" -lt "$(echo $edt | tr -d '_')00" ]; then 
        echo $filepath >> $configdir/$tfile
        echo $(basename $filepath)
    fi
    if [ "$filedt" -ge "$(echo $edt | tr -d '_')00" ]; then
        break
    fi
done

mv $configdir/$tfile $configdir/$inputfile

tfile=tmp_${stt}_${edt}.txt

if [[ "$(ls $tempdir/* | sort | head -n 1 | xargs basename)" == "wrfout"* ]]; then

    echo
    echo Selecting temperature files for analysis in range $stt to $edt...
    echo
    sleep 3 

    for filepath in $(ls $tempdir/* | sort); do
         
        file=$(basename $filepath)
        filedt=$(echo $file | cut -d '_' -f4 | tr -d '-')$(echo $file | cut -d '_' -f5 | cut -d '.' -f1 | tr -d ':')
        if [ "$filedt" -ge "$(echo $stt | tr -d '_')00" ] && [ "$filedt" -le "$(echo $edt | tr -d '_')00" ]; then 
            echo $filepath >> $configdir/$tfile
            echo $(basename $filepath)
        fi
        if [ "$filedt" -gt "$(echo $edt | tr -d '_')00" ]; then
            break
        fi

    done

    mp=$(head -n 1 $configdir/$tfile | xargs basename | cut -d '_' -f2)
    for ((ii=0;ii<${#mpopts[@]};ii++)); do
       [[ "${mpopts[ii]}" = "$mp" ]] && break
    done
    mpname=$(echo "${mpnames[ii]}")
    tempsrc=$mpname
    tempfile=temp_${tempsrc}_${stt}_${edt}.txt
    snd_on='False'
    wrft_on='True'

else

    echo
    echo Selecting sounding files at closest times to range $stt to $edt...
    echo
    sleep 3

    declare -a beffiles=()
    declare -a infiles=()
    declare -a aftfiles=()

    for filepath in $(ls $tempdir/* | sort); do

        file=$(basename $filepath)
        filedt=$(echo $file | cut -d '_' -f2 | tr -d '-')$(echo $file | cut -d '_' -f3 | cut -d '.' -f1 | tr -d ':')
       
        if [ "$filedt" -lt "$(echo $stt | tr -d '_')00" ]; then
            beffiles+=($filepath)
        elif [ "$filedt" -ge "$(echo $stt | tr -d '_')00" ] && [ "$filedt" -le "$(echo $edt | tr -d '_')00" ]; then
            infiles+=($filepath)
        elif [ "$filedt" -gt "$(echo $edt | tr -d '_')00" ]; then 
            aftfiles+=($filepath)
        fi
        
    done

    echo $(basename ${beffiles[${#beffiles[@]}-1]})
    echo ${beffiles[${#beffiles[@]}-1]} > $configdir/$tfile
    for filep in ${infiles[@]}; do
        echo $(basename $filep)
        echo $filep >> $configdir/$tfile
    done
    echo $(basename ${aftfiles[0]})
    echo ${aftfiles[0]} >> $configdir/$tfile
    
    tempsrc='uwyo-'$(basename $tempdir | tr '[:upper:]' '[:lower:]')
    tempfile=temp_${tempsrc}_${stt}_${edt}.txt
    snd_on='True'
    wrft_on='False'

fi

mv $configdir/$tfile $configdir/$tempfile

if [ ! -z $doppdir ]; then

    echo
    echo Selecting dual-Doppler files at closest times to range $stt to $edt...
    echo
    sleep 3

    dd_on='True'
    doppfile=dopp_${data}_${stt}_${edt}.txt
    
    tfile=dopp.txt

    for filepath in $(ls $doppdir/* | sort); do
        file=$(basename $filepath) 
        if [[ "$(ls $raddir/* | head -n 1 | xargs basename)" == "wrfout"* ]]; then
            filedt=$(echo $file | cut -d '_' -f4 | tr -d '-')$(echo $file | cut -d '_' -f5 | cut -d '.' -f1 | tr -d ':')
        else
            #filedt=$(echo $file | cut -d '_' -f2)$(echo $file | cut -d '_' -f3)
            filedt=$(echo $file | cut -d '_' -f3)$(echo $file | cut -d '_' -f4)
        fi
        if [ "$filedt" -ge "$(echo $stt | tr -d '_')00" ] && [ "$filedt" -lt "$(echo $edt | tr -d '_')00" ]; then 
            echo $filepath >> $configdir/$tfile
            echo $(basename $filepath)
        fi
        if [ "$filedt" -ge "$(echo $edt | tr -d '_')00" ]; then
            break
        fi
    done

    mv $configdir/$tfile $configdir/$doppfile

else
    dd_on='False'
fi

if [ -z $simdirs ]; then
    if [[ "$data" == "obs" ]]; then
        fold="obs"
    else
        fold=$mpname
    fi

else

    declare -a inputfiles2=()
    declare -a configfiles2=()
    declare -a allmps=()

    for simdir in ${simdirs[@]}; do
        
        echo
        echo Selecting wrfout files for analysis in range $stt to $edt...

        fold='obsvwrf'
        mp2=$(ls $simdir/* | head -n 1 | xargs basename | cut -d '_' -f2)
        for ((ii=0;ii<${#mpopts[@]};ii++)); do
           [[ "${mpopts[ii]}" = "$mp2" ]] && break
        done
        mpname2=$(echo ${mpnames[ii]})
        echo MP=$(echo $mpname2 | tr '[:lower:]' '[:upper:]')
        echo

        inputfile2=input_${mpname2}_${stt}_${edt}.txt
        configfile2=config_${mpname2}_${stt}_${edt}.txt

        inputfiles2+=($configdir/$inputfile2)
        configfiles2+=($configdir/$configfile2)
        allmps+=($mpname2)

        tfile=input_${stt}_${edt}.txt

        for filepath in $(ls $simdir/* | sort); do
            file=$(basename $filepath)
            filedt=$(echo $file | cut -d '_' -f4 | tr -d '-')$(echo $file | cut -d '_' -f5 | cut -d '.' -f1 | tr -d ':')
            if [ "$filedt" -ge "$(echo $stt | tr -d '_')00" ] && [ "$filedt" -lt "$(echo $edt | tr -d '_')00" ]; then 
                echo $filepath >> $configdir/$tfile
                echo $(basename $filepath)
            fi
            if [ "$filedt" -ge "$(echo $edt | tr -d '_')00" ]; then
                break
            fi
        done

        mv $configdir/$tfile $configdir/$inputfile2

    done

fi

outfigdir=outputfig/${fold}_temp${tempsrc}_${station}_${stt}_${edt}
outrrdir=$(cd $raddir/../../ && pwd)/radar_rainrates/$station
mkdir -p $outfigdir $outrrdir

if [ -z $simdir ]; then

    echo
    echo Creating config file for iPOLARRIS...
    echo
    sleep 3

    template=$ipoldir/${agency}_${data}_config.txt
    cp $template $configdir/$configfile

    sed -i '' "s/^type ==.*/type == $data == # Type of input data: 'obs' OR 'wrf' (obs + simulated)/g" $configdir/$configfile
    if [[ "$data" == "obs" ]]; then
        sed -i '' "s/.*mphys ==.*/mphys == $data == # Type of microphysics used in model: 'obs' OR '<scheme>' if type = 'wrf'/g" $configdir/$configfile
    else
        sed -i '' "s/.*mphys ==.*/mphys == $mpname == # Type of microphysics used in model: 'obs' OR '<scheme>' if type = 'wrf'/g" $configdir/$configfile
    fi
    sed -i '' "s/.*ptype ==.*/ptype == '$ptype' == # Output figure file extenstion (i.e. png, jpg, mp4, ...)/g" $configdir/$configfile
    sed -i '' "s/.*sdatetime ==.*/sdatetime == '$(echo $stt | tr '_' '-')' == # Start time of analysis of interest/g" $configdir/$configfile
    sed -i '' "s/.*edatetime ==.*/edatetime == '$(echo $edt | tr '_' '-')' == # End time of analysis of interest/g" $configdir/$configfile
    sed -i '' "s%.*rfiles ==.*%rfiles == '$configdir/$inputfile' == # Path to list of radar files to read in%g" $configdir/$configfile
    if [[ "$wrft_on" == "True" ]]; then
        sed -i '' "s%.*wfiles ==.*%wfiles == '$configdir/$tempfile' == # Path to list of WRF temperature files to read in%g" $configdir/$configfile
    elif [[ "$snd_on" == "True" ]]; then
        sed -i '' "s%.*sfiles ==.*%sfiles == '$configdir/$tempfile' == # Path to list of sounding files to read in%g" $configdir/$configfile
    fi
    if [[ "$dd_on" == "True" ]]; then
        sed -i '' "s%.*dfiles ==.*%dfiles == '$configdir/$doppfile' == # Path to list of dual-Doppler files to read in%g" $configdir/$configfile
    fi
    if [[ "$data" == "obs" ]]; then
        sed -i '' "s/.*exper ==.*/exper == $station == # Radar location/g" $configdir/$configfile
    else
        sed -i '' "s/.*exper ==.*/exper == $station-$(echo $mpname | tr '[:lower:]' '[:upper:]') == # Radar location/g" $configdir/$configfile
    fi
    sed -i '' "s/lat ==  == #.*/lat == $latcen == # Latitude of the radar station/g" $configdir/$configfile
    sed -i '' "s/lon ==  == #.*/lon == $loncen == # Longitude of the radar station/g" $configdir/$configfile
    sed -i '' "s/costhresh ==.*/costhresh == $maxelev == # Elevation angle threshold for cone of silence/g" $configdir/$configfile
    sed -i '' "s%.*image_dir ==.*%image_dir == '$outfigdir/' == # Output figure directory%g" $configdir/$configfile
    sed -i '' "s%.*rr_dir ==.*%rr_dir == '$outrrdir/' == # Output rain rate netcdf directory%g" $configdir/$configfile
    sed -i '' "s/.*dd_on ==.*/dd_on == $dd_on == # Doppler gridded velocity on/g" $configdir/$configfile
    sed -i '' "s/.*snd_on ==.*/snd_on == $snd_on == # Sounding temperature on/g" $configdir/$configfile
    sed -i '' "s/.*wrft_on ==.*/wrft_on == $wrft_on == # WRF temperature on/g" $configdir/$configfile

    echo Running iPOLARRIS...
    sleep 3

    python run_ipolarris.py $configdir/$configfile

else

    echo
    echo Creating OBS and SIM config files for iPOLARRIS...
    echo
    sleep 3

    echo OBS
    template=$ipoldir/${agency}_obs_config.txt
    cp $template $configdir/$configfile

    sed -i '' "s/^type ==.*/type == obs == # Type of input data: 'obs' OR 'wrf' (obs + simulated)/g" $configdir/$configfile
    sed -i '' "s/.*mphys ==.*/mphys == obs == # Type of microphysics used in model: 'obs' OR '<scheme>' if type = 'wrf'/g" $configdir/$configfile
    sed -i '' "s/.*ptype ==.*/ptype == '$ptype' == # Output figure file extenstion (i.e. png, jpg, mp4, ...)/g" $configdir/$configfile
    sed -i '' "s/.*sdatetime ==.*/sdatetime == '$(echo $stt | tr '_' '-')' == # Start time of analysis of interest/g" $configdir/$configfile
    sed -i '' "s/.*edatetime ==.*/edatetime == '$(echo $edt | tr '_' '-')' == # End time of analysis of interest/g" $configdir/$configfile
    sed -i '' "s%.*rfiles ==.*%rfiles == '$configdir/$inputfile' == # Path to list of radar files to read in%g" $configdir/$configfile
    if [[ "$wrft_on" == "True" ]]; then
        sed -i '' "s%.*wfiles ==.*%wfiles == '$configdir/$tempfile' == # Path to list of WRF temperature files to read in%g" $configdir/$configfile
    elif [[ "$snd_on" == "True" ]]; then
        sed -i '' "s%.*sfiles ==.*%sfiles == '$configdir/$tempfile' == # Path to list of sounding files to read in%g" $configdir/$configfile
    fi
    if [[ "$dd_on" == "True" ]]; then
        sed -i '' "s%.*dfiles ==.*%dfiles == '$configdir/$doppfile' == # Path to list of dual-Doppler files to read in%g" $configdir/$configfile
    fi
    sed -i '' "s/.*exper ==.*/exper == $station == # Radar location/g" $configdir/$configfile
    sed -i '' "s/lat ==  == #.*/lat == $latcen == # Latitude of the radar station/g" $configdir/$configfile
    sed -i '' "s/lon ==  == #.*/lon == $loncen == # Longitude of the radar station/g" $configdir/$configfile
    sed -i '' "s%.*image_dir ==.*%image_dir == '$outfigdir/' == # Output figure directory%g" $configdir/$configfile
    sed -i '' "s%.*rr_dir ==.*%rr_dir == '$outrrdir/' == # Output rain rate netcdf directory%g" $configdir/$configfile
    sed -i '' "s/.*dd_on ==.*/dd_on == $dd_on == # Doppler gridded velocity on/g" $configdir/$configfile
    sed -i '' "s/.*snd_on ==.*/snd_on == $snd_on == # Sounding temperature on/g" $configdir/$configfile
    sed -i '' "s/.*wrft_on ==.*/wrft_on == $wrft_on == # WRF temperature on/g" $configdir/$configfile

    template2=$ipoldir/${agency}_wrf_config.txt
    
    for ((ii=0;ii<${#allmps[@]};ii++)); do
        
        echo MP=$(echo ${allmps[ii]} | tr '[:lower:]' '[:upper:]')
        cp $template2 ${configfiles2[ii]}

        sed -i '' "s/^type ==.*/type == wrf == # Type of input data: 'obs' OR 'wrf' (obs + simulated)/g" ${configfiles2[ii]}
        sed -i '' "s/.*mphys ==.*/mphys == $(echo ${allmps[ii]}) == # Type of microphysics used in model: 'obs' OR '<scheme>' if type = 'wrf'/g" ${configfiles2[ii]}
        sed -i '' "s/.*ptype ==.*/ptype == '$ptype' == # Output figure file extenstion (i.e. png, jpg, mp4, ...)/g" ${configfiles2[ii]}
        sed -i '' "s/.*sdatetime ==.*/sdatetime == '$(echo $stt | tr '_' '-')' == # Start time of analysis of interest/g" ${configfiles2[ii]}
        sed -i '' "s/.*edatetime ==.*/edatetime == '$(echo $edt | tr '_' '-')' == # End time of analysis of interest/g" ${configfiles2[ii]}
        sed -i '' "s%.*rfiles ==.*%rfiles == '${inputfiles2[ii]}' == # Path to list of radar files to read in%g" ${configfiles2[ii]}
        sed -i '' "s%.*wfiles ==.*%wfiles == '${inputfiles2[ii]}' == # Path to list of WRF temperature files to read in%g" ${configfiles2[ii]}
        if [[ "$dd_on" == "True" ]]; then
            sed -i '' "s%.*dfiles ==.*%dfiles == '${inputfiles2[ii]}' == # Path to list of dual-Doppler files to read in%g" ${configfiles2[ii]}
        fi
        sed -i '' "s/.*exper ==.*/exper == $station-$(echo ${allmps[ii]} | tr '[:lower:]' '[:upper:]') == # Radar location/g" ${configfiles2[ii]}
        sed -i '' "s/lat ==  == #.*/lat == $latcen == # Latitude of the radar station/g" ${configfiles2[ii]}
        sed -i '' "s/lon ==  == #.*/lon == $loncen == # Longitude of the radar station/g" ${configfiles2[ii]}
        sed -i '' "s/costhresh ==.*/costhresh == $maxelev == # Elevation angle threshold for cone of silence/g" ${configfiles2[ii]}
        sed -i '' "s%.*image_dir ==.*%image_dir == '$outfigdir/' == # Output figure directory%g" ${configfiles2[ii]}
        sed -i '' "s%.*rr_dir ==.*%rr_dir == '$outrrdir/' == # Output rain rate netcdf directory%g" ${configfiles2[ii]}
        sed -i '' "s/.*dd_on ==.*/dd_on == $dd_on == # Doppler gridded velocity on/g" ${configfiles2[ii]}
        sed -i '' "s/.*snd_on ==.*/snd_on == False == # Sounding temperature on/g" ${configfiles2[ii]}
        sed -i '' "s/.*wrft_on ==.*/wrft_on == True == # WRF temperature on/g" ${configfiles2[ii]}

    done

    echo
    echo Running iPOLARRIS...
    sleep 3

    python run_ipolarris.py $configdir/$configfile $(printf "%s " ${configfiles2[@]})

fi
