# Astronomy math


# Build
	Open solution astronomy_sdk_vs2015.sln and build

# Get data
	cd to ./bin
	run get_data.cmd
	run get_x64.cmd if build x64
	run get_x86.cmd if build x86
	
	
	
	
save
//
stages:
  - build

job:
  stage: build
  script:
  - git submodule update --init
  - ''
  - echo "Release build..."
  - ''
  - 'Start-Process "cmd.exe" "/c build_sdk.bat"'
  - ''
  - echo "Build success..."
  tags: 
  - my-tag
  except:
  - tags