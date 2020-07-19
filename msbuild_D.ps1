function buildVS
{
    param
    (
        [parameter(Mandatory=$true)]
        [String] $path,

        [parameter(Mandatory=$false)]
        [bool] $nuget = $true,
        
        [parameter(Mandatory=$false)]
        [bool] $clean = $true
    )
    process
    {
        $msBuildExe = 'C:\Program Files (x86)\MSBuild\14.0\Bin\msbuild.exe'

        if ($nuget) {
            Write-Host "Restoring NuGet packages" -foregroundcolor green
            nuget restore "$($path)"
        }

        if ($clean) {
            Write-Host "Cleaning $($path)" -foregroundcolor green
            & "$($msBuildExe)" "$($path)" /t:Clean /m
        }

        Write-Host "Building $($path)" -foregroundcolor green
        & "$($msBuildExe)" "$($path)" /t:ReBuild /property:Configuration=Debug /property:Platform=x64 /m
		$MyLastExitCode = $LastExitCode
		
        Write-Host ("BUILD END")
        if( $MyLastExitCode -eq 0 )	
        {	
			Write-Host ("BUILD OK")
        }
        else
        {			
			Write-Host ("BUILD ERROR")	
            Write-Error ( $MyLastExitCode )
			exit 1
        }
    }
}

chcp 850
buildVS -path .\astronomy_sdk_vs2015.sln -nuget $false -clean $true