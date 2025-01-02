
cd external
wget https://elsevier.digitalcommonsdata.com/public-files/datasets/yc9vv7rwyj/files/be55b15b-6d9c-4f0e-b5f7-09b300cfb806/file_downloaded
tar -xzf file_downloaded
rm file_downloaded

echo "Building Zeal..."
cd Zeal
make -j
echo "Finished building Zeal..."
echo ""
cd ../../KIM/

mkdir build
cd build
cmake ..
make -j
