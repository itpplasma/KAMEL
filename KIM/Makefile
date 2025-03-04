.PHONY: all ninja test install clean
all: ninja

build/build.ninja:
	wget https://elsevier.digitalcommonsdata.com/public-files/datasets/yc9vv7rwyj/files/be55b15b-6d9c-4f0e-b5f7-09b300cfb806/file_downloaded -O Zeal.tar.gz
	tar -xzf Zeal.tar.gz
	rm Zeal.tar.gz
	mv Zeal ../external/
	cmake -GNinja -S . -B build

ninja: build/build.ninja
	cmake --build build

clean:
	rm -rf ../external/Zeal
	rm -rf build