call "%VS140COMNTOOLS%\vsvars32.bat"

devenv "astronomy_sdk_vs2015.sln" /rebuild "Release|x64"

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