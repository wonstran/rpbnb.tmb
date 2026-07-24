from io import BytesIO
from pathlib import Path
import sys

from pypdf import PdfReader, PdfWriter
from reportlab.pdfgen import canvas


def stamp(input_path: Path, output_path: Path) -> None:
    reader = PdfReader(str(input_path))
    writer = PdfWriter()
    total = len(reader.pages)

    for number, page in enumerate(reader.pages, start=1):
        width = float(page.mediabox.width)
        height = float(page.mediabox.height)
        packet = BytesIO()
        overlay_canvas = canvas.Canvas(packet, pagesize=(width, height))
        overlay_canvas.setFont("Helvetica", 8)
        overlay_canvas.drawCentredString(
            width / 2, 18, f"rpbnb.tmb Reference Manual - {number} of {total}"
        )
        overlay_canvas.save()
        packet.seek(0)
        page.merge_page(PdfReader(packet).pages[0])
        writer.add_page(page)

    with output_path.open("wb") as handle:
        writer.write(handle)


if __name__ == "__main__":
    if len(sys.argv) != 3:
        raise SystemExit("usage: add_pdf_page_numbers.py INPUT.pdf OUTPUT.pdf")
    stamp(Path(sys.argv[1]), Path(sys.argv[2]))
