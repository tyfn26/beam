using GLMakie
using GeometryBasics
using Colors

# --- 1. Ana Pencere ve Tema Ayarları ---
# Ansys Mechanical'ın koyu temasını ve fontlarını taklit ediyoruz
set_theme!(theme_dark(), font = "Segoe UI")

function launch_ansys_clone()
    # Pencere boyutu ve başlığı
    fig = Figure(size = (1400, 900), backgroundcolor = RGBf(0.2, 0.2, 0.2))
    
    # --- 2. Yerleşim Düzeni (Layout Grid) ---
    # Sol: Proje Ağacı (Sidebar) | Orta: 3D Viewport | Sağ: Renk Skalası
    # Alt: Mesaj/Konsol
    
    # Grid tanımları
    sidebar_grid = fig[1, 1] = GridLayout()
    viewport_grid = fig[1, 2]
    colorbar_grid = fig[1, 3]
    msg_grid = fig[2, :]

    # Sütun genişliklerini ayarla (Sidebar sabit, Viewport esnek)
    colsize!(fig.layout, 1, Fixed(250)) 
    colsize!(fig.layout, 3, Fixed(100))
    rowsize!(fig.layout, 2, Fixed(80))

    # --- 3. Sol Panel: Proje Ağacı (Simulation Tree) ---
    # Arka plan kutusu
    Box(sidebar_grid[1:10, 1], color = RGBf(0.9, 0.9, 0.9))
    
    # Başlık
    Label(sidebar_grid[1, 1], "Project Outline", textsize = 20, color = :black, font = :bold, halign = :left, padding = (10,0,10,0))
    
    # Ağaç Elemanları (Buton görünümlü etiketler)
    items = ["Model (A4)", "  > Geometry", "  > Materials", "  > Mesh", "Static Structural (A5)", "  > Analysis Settings", "  > Loads & Supports", "Solution (A6)", "  > Total Deformation", "  > Equivalent Stress"]
    
    # Basit bir döngü ile ağacı oluşturuyoruz
    for (i, item) in enumerate(items)
        color = startswith(item, "  ") ? :black : :blue # Ana başlıklar mavi
        font_w = startswith(item, "  ") ? :regular : :bold
        Label(sidebar_grid[i+1, 1], item, color = color, halign = :left, textsize = 16, font = font_w, padding = (15,0,0,0))
    end
    
    # Boşluğu doldurmak için görünmez eleman
    Label(sidebar_grid[11, 1], " ", tellheight=true)

    # --- 4. Detaylar Paneli (Sol Alt) ---
    Label(sidebar_grid[12, 1], "Details of 'Stress'", textsize = 18, color = :black, font = :bold, halign = :left, padding = (10,10,0,0))
    details_text = "Type: Equivalent (von-Mises)\nUnit: MPa\nMax: 250.4 MPa\nMin: 0.02 MPa\nElement Order: Quadratic"
    Label(sidebar_grid[13, 1], details_text, color = :black, halign = :left, textsize = 14, padding = (15,0,20,0))

    # --- 5. Orta Panel: 3D Görüntüleme (Viewport) ---
    # LScene, 3D etkileşim (döndürme/zoom) sağlayan asıl bileşendir.
    # Ansys mavisi arka plan
    ax3d = LScene(viewport_grid, scenekw = (backgroundcolor = RGBf(0.2, 0.2, 0.35), clear=true))
    
    # --- FEA Modeli Oluşturma (Bir Kiriş/Beam Simülasyonu) ---
    # Basit bir dikdörtgen prizma (Beam) oluşturuyoruz
    x = range(0, 10, length=30)
    y = range(0, 2, length=10)
    z = range(0, 2, length=10)
    
    # Stres verisi uyduralım (Ankastre kiriş gibi: kökte stres çok, uçta az)
    # meshscatter veya surface yerine 'mesh' ve 'color' kullanacağız.
    
    # Noktaları ve renk değerlerini tanımla
    points = [Point3f(xi, yi, zi) for xi in x, yi in y, zi in z]
    stress_values = [250 * (1 - xi/10)^2 + rand()*10 for xi in x, yi in y, zi in z] # Kök tarafta yüksek stres
    
    # Mesh çizimi (Volume/Contour Plot benzeri)
    # Not: Gerçek bir FEA mesh'i için 'mesh' fonksiyonu kullanılır, burada görsel şölen için hacimsel scatter kullanıyoruz.
    plt = meshscatter!(ax3d, vec(points), 
        color = vec(stress_values), 
        colormap = :turbo, # Ansys tarzı gökkuşağı renkleri
        markersize = 0.15,
        shading = FastShading # Işıklandırma ekle
    )
    
    # Eksenleri gizle (Ansys'de genelde model boştadır)
    # hidedecorations!(ax3d) # İstersen açabilirsin

    # --- 6. Sağ Panel: Renk Skalası (Legend) ---
    Colorbar(colorbar_grid[1,1], plt, label = "Equivalent Stress (MPa)", height = Relative(0.6), width = 20, ticklabelsize = 14)

    # --- 7. Alt Panel: Mesaj/Output ---
    Box(msg_grid, color = RGBf(0.95, 0.95, 0.95))
    msg_txt = "Messages:\n> Solver started...\n> Global matrix assembled.\n> Solution converged in 4 iterations.\n> Result file loaded successfully."
    Label(msg_grid[1, 1], msg_txt, color = :black, font = "Consolas", textsize = 12, halign = :left, valign = :top, padding = (10,10,10,10))

    # Pencereyi göster
    display(fig)
end

# Programı çalıştır
launch_ansys_clone()
