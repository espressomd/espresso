#!/bin/sh

echo "GENERATING AND OPTIMIZING PDFs..."
inkscape -z --export-pdf=logo-noopt.pdf logo.svg
pdfopt logo-noopt.pdf logo.pdf
inkscape -z --export-pdf=logo-small-noopt.pdf logo-small.svg
pdfopt logo-small-noopt.pdf logo-small.pdf

rm -f logo*-noopt.pdf
ls -l logo.pdf logo-small.pdf

echo
echo "GENERATING 500x500 PNG..."
inkscape -z --export-png=logo_500x500.png --export-width=500 --export-height=500 logo.svg
echo "GENERATING 100x100 PNG..."
inkscape -z --export-png=logo_100x100.png --export-width=100 --export-height=100 logo.svg
echo "GENERATING 32x32 PNG..."
inkscape -z --export-png=logo_32x32.png --export-width=32 --export-height=32 logo-small.svg
echo "GENERATING 48x48 PNG..."
inkscape -z --export-png=logo_48x48.png --export-width=48 --export-height=48 logo-small.svg
echo

ls -l \
  logo_500x500.png \
  logo_100x100.png \
  logo_48x48.png \
  logo_32x32.png

