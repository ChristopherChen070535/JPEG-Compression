INCLUDE = -I.

all: build_encoder build_decoder

build_encoder:
	gcc -O3 -o encoder_411086013.exe encoder_411086013.c -lm

build_decoder:
	gcc -O3 -o decoder_411086013.exe decoder_411086013.c -lm

demo0:
	./encoder_411086013.exe 0 KimberlyNCat.bmp R.txt G.txt B.txt dim.txt
	./decoder_411086013.exe 0 ResKimberlyNCat.bmp R.txt G.txt B.txt dim.txt

demo1:
	./encoder_411086013.exe 1 KimberlyNCat.bmp Qt_Y.txt Qt_Cb.txt Qt_Cr.txt dim.txt qF_Y.raw qF_Cb.raw qF_Cr.raw eF_Y.raw eF_Cb.raw eF_Cr.raw
	./decoder_411086013.exe 1 QResKimberlyNCat.bmp KimberlyNCat.bmp Qt_Y.txt Qt_Cb.txt Qt_Cr.txt dim.txt qF_Y.raw qF_Cb.raw qF_Cr.raw
	./decoder_411086013.exe 1 ResKimberlyNCat1.bmp Qt_Y.txt Qt_Cb.txt Qt_Cr.txt dim.txt qF_Y.raw qF_Cb.raw qF_Cr.raw eF_Y.raw eF_Cb.raw eF_Cr.raw

demo2:
	./encoder_411086013.exe 2 KimberlyNCat.bmp ascii rle_code.txt
	./decoder_411086013.exe 2 QResKimberlyNCat2.bmp ascii rle_code.txt

clean:
	rm -f encoder_411086013.exe decoder_411086013.exe R.txt G.txt B.txt dim.txt Qt_Y.txt Qt_Cb.txt Qt_Cr.txt qF_Y.raw qF_Cb.raw qF_Cr.raw eF_Y.raw eF_Cb.raw eF_Cr.raw QResKimberlyNCat.bmp ResKimberlyNCat.bmp ResKimberlyNCat1.bmp QResKimberlyNCat2.bmp rle_code.txt

.PHONY: build_encoder build_decoder demo0 demo1 clean
