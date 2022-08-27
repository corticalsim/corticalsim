# modelCMT (cortical microtubule dynamics)

![alt text](https://github.com/HiBandan/modelCMT/blob/main/modules/logo/logomodelCMT-1.0.png)


This is a program for simulating cortical microtubule dynamics (CMT) on experimentally extracted microscopic images of cell 

Application: 

https://www.sciencedirect.com/science/article/pii/S0960982218309242

https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1005959

## Note

Testing and contributing is very welcome, especially if you can contribute with new algorithms and features.

## WINDOWS

  ### install gcc/g++ 
  
  download Mingw: https://nuwen.net/mingw.html
  
      run mingw-xx.exe -> this will create a folder MinGw
      
      copy MinGw to C:\
       
      add environment variable: C:\MinGW or C:\MinGW\bin

  ### install eigen:
  
  download eigen: https://drive.google.com/file/d/1VFLlKJI9EaSP9VUaeV8vmPQYkqdIIJPg/view?usp=sharing
  
      extract "eigen-3.4.0" folder and rename it to "Eigen"
  
      copy "Eigen" folder to C:\Program Files 
  
  ### install boost:
  
  download boost: https://drive.google.com/file/d/1apu5_am2kJj7HvXJNPhn30k_ryi3gJDf/view?usp=sharing

      extract "boost_1_77_0" folder and rename it to "Boost"
  
      copy "Boost" to C:\Program Files 
  
  ### get the source code
  
  Option 1: Download the package from the Github: https://github.com/HiBandan/modelCMT/archive/refs/heads/main.zip


  Option 2: Clone from terminal: 
  
    git clone https://github.com/HiBandan/modelCMT.git

  ### install Code::Blocks: 
  https://code-blocks.en.uptodown.com/windows
    
  ### run 
  
    open (with Code::Blocks): modelCMT.cbp
  
    Projects -> modelCMT.cbp -> add files
        -> cppFiles (all files except 'inline.cpp')
        -> headFiles (all files)
        -> paramFiles (all files)
        
    settings -> compiler -> search directories (compiler) -> add 
  
        -> C:\Program Files\Eigen 
        -> C:\Program Files\Boost 

    project -> 
    
      set program’s arguments: ./paramFiles/parameters_ARRAY.txt
      build options -> search directories -> compiler -> add -> cppFiles 
      build options -> search directories -> compiler -> add -> headFiles
  
    run 

## LINUX

  ### install gcc/g++ 
  
  download Mingw: https://nuwen.net/mingw.html
  
       sudo apt install build-essential

  ### install eigen:
  
  download eigen: https://drive.google.com/file/d/1VFLlKJI9EaSP9VUaeV8vmPQYkqdIIJPg/view?usp=sharing
  
      extract "eigen-3.4.0" folder and rename it to "Eigen"
  
      copy "Eigen" folder to /usr/local/include
  
  ### install boost:
  
  download boost: https://drive.google.com/file/d/1apu5_am2kJj7HvXJNPhn30k_ryi3gJDf/view?usp=sharing

      extract "boost_1_77_0" folder and rename it to "Boost"
  
      copy "Boost" to /usr/local/include

  ### get the source code
  
  Option 1: Download the package from the Github: https://github.com/HiBandan/modelCMT/archive/refs/heads/main.zip

  Option 2: Clone from terminal: 
  
    git clone https://github.com/HiBandan/modelCMT.git
    
  ### install Code::Blocks: 
  
  sudo add-apt-repository ppa:damien-moore/codeblocks-stable

  sudo apt update

  sudo apt install codeblocks codeblocks-contrib
 
  ### run 

    open (with code-block): modelCMT.cbp
    
    Projects -> modelCMT.cbp -> add files
      -> cppFiles (all files except 'inline.cpp')
      -> headFiles (all files)
      -> paramFiles (all files)
  
    settings -> compiler -> search directories (compiler) -> add 
  
        -> /usr/local/include/Eigen
        -> /usr/local/include/Boost 
       
    project -> 
    
      set program’s arguments: ./paramFiles/parameters_ARRAY.txt
      build options -> search directories -> compiler -> add -> cppFiles 
      build options -> search directories -> compiler -> add -> headFiles
  
    run 
