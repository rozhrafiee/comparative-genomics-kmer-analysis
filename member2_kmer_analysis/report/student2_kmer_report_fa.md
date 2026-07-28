<div dir="rtl" lang="fa">

# ۱. عنوان پروژه

**ژنومیک مقایسه‌ای: تحلیل K-mer و محتوای GC در پنج ارگانیسم — بخش دانشجوی ۲**

این گزارش فقط بخش تحلیل K-mer، Shannon Entropy، ماتریس شباهت و فاصله، PCA، خوشه‌بندی سلسله‌مراتبی، دندروگرام، جستجوی TATA box، تحلیل آماری (همبستگی، ANOVA، Kruskal–Wallis، Mann–Whitney) و نمودارهای مربوط به همین تحلیل‌ها را پوشش می‌دهد.

---

## ۲. چکیده

پنج اسمبلی ژنومی NCBI (`ncbi_dataset1` تا `ncbi_dataset5`) با پایپ‌لاین قابل بازتولید در پوشه `member2_kmer_analysis` تحلیل شدند. هویت گونه‌ها از `data_summary.tsv` استخراج شد: *Caenorhabditis elegans*، *Neurospora crassa*، *Yarrowia lipolytica*، *Dictyostelium discoideum* و *Eremothecium coryli*. فهرست گونه‌های Brief فقط نمونه بود؛ تحلیل روی همین پنج اسمبلی واقعی انجام شد.

برای هر گونه، همه کروموزوم‌ها/اسکفولدها در یک پروفایل k-mer تجمیع شدند؛ تحلیل اصلی با k=۴ و جدول پرتکرارترین‌ها برای k=۳، ۴ و ۵ انجام شد. هر پنج ژنوم همه ۲۵۶ نوع ۴-mer ممکن را داشتند و هیچ k-mer گونه‌اختصاصی در k=۳–۵ یافت نشد. بیشترین Shannon entropy مربوط به *N. crassa* (۷٫۹۶۵ بیت) و بیشترین `kmer_diversity` مربوط به *E. coryli* بود. نزدیک‌ترین جفت از نظر فاصله اقلیدسی *N. crassa* و *Y. lipolytica* (۰٫۰۱۵۴) و دورترین جفت *Y. lipolytica*–*D. discoideum* (۰٫۱۰۲۳) بود. در PCA، PC1 حدود ۷۰٫۴٪ واریانس را توضیح داد.

در بخش آماری، همبستگی سطح‌گونه (n=۵) بین اندازه ژنوم و entropy ضعیف و غیرمعنادار بود؛ در مقابل، همبستگی GC% (خوانده‌شده فقط‌خواندنی از خروجی دانشجوی ۱) با Shannon entropy قوی و معنادار گزارش شد (Pearson r≈۰٫۹۶۶، p≈۰٫۰۰۷) ولی با توان آماری پایین به‌خاطر n کوچک. آزمون‌های گروهی روی رکوردهای جدا (کروموزوم/اسکفولد) نشان داد entropy و k-mer diversity بین گونه‌ها تفاوت معنادار دارند (ANOVA و Kruskal–Wallis با p≪۰٫۰۰۱)، با احتیاط شبه-تکرار بودن رکوردهای داخل یک ژنوم.

جستجوی TATA box در ۵۰۰ bp اول توالی تجمیع‌شده، موتیف `TATAAA` را در *N. crassa* (۵ محل) و *D. discoideum* (۲ محل) یافت. تحلیل GC و اندازه ژنوم متعلق به بخش دانشجوی ۱ است و در این ماژول فقط برای همبستگی اختیاری خوانده می‌شود.

---

## ۳. مقدمه

K-mer یک زیررشته با طول ثابت k از توالی DNA است. توزیع فراوانی k-merها یک «اثر انگشت» ترکیبی از ژنوم می‌سازد که بدون هم‌ترازسازی کامل برای مقایسه گونه‌ها، خوشه‌بندی و کاهش بُعد استفاده می‌شود. Shannon entropy یکنواختی این توزیع را می‌سنجد. فاصله اقلیدسی و cosine similarity شباهت پروفایل‌ها را کمی می‌کنند. PCA و خوشه‌بندی سلسله‌مراتبی ساختار روابط را نشان می‌دهند. آزمون‌های آماری (Pearson/Spearman، ANOVA، Kruskal–Wallis، Mann–Whitney) روابط و تفاوت‌های گروهی را بررسی می‌کنند. TATA box یک موتیف پروموتری یوکاریوتی است که در این بخش به‌صورت اسکن دنباله در ابتدای توالی بررسی شده است.

---

## ۴. اهداف بخش دوم پروژه

- شمارش و نرمال‌سازی k-merها برای هر ژنوم؛
- محاسبه Shannon entropy توزیع k-mer؛
- ساخت ماتریس فاصله اقلیدسی و cosine similarity؛
- اجرای PCA و خوشه‌بندی سلسله‌مراتبی (Ward) همراه با دندروگرام؛
- جستجوی موتیف TATA box در ناحیه ابتدایی توالی؛
- تحلیل آماری: همبستگی سطح‌گونه و آزمون‌های گروهی روی رکوردها؛
- تولید heatmap، PCA، dendrogram، entropy، و نمودارهای آمار؛
- آماده‌سازی گزارش فارسی و راهنمای ارائه.

---

## ۵. مواد و روش‌ها

اسکریپت‌های پوشه `scripts/` با Python، Biopython، NumPy، Pandas، SciPy، Scikit-learn، Matplotlib و Seaborn اجرا شدند. ورودی پیش‌فرض پوشه نسبی `../data` است. دیتاست‌های `ncbi_dataset*` به‌صورت خودکار کشف می‌شوند؛ برای هر پوشه، فایل RefSeq (`GCF_*_genomic.fna`) در صورت وجود و در غیر این صورت GenBank (`GCA_*_genomic.fna`) انتخاب می‌شود.

توالی‌ها به حروف بزرگ تبدیل و فقط بازهای A/C/G/T در شمارش k-mer نگه داشته می‌شوند. فراوانی هر k-mer بر مجموع k-merهای معتبر نرمال می‌شود تا پروفایل مستقل از طول ژنوم باشد.

برای آزمون‌های گروهی (ANOVA / Kruskal–Wallis / Mann–Whitney)، تابع `load_sequences_by_record` هر رکورد FASTA را جدا نگه می‌دارد تا برای هر گونه چند نمونه (کروموزوم/اسکفولد) وجود داشته باشد. خروجی‌های جدولی (شامل آمار) در `results/` و شکل‌ها در `results/figures/` ذخیره می‌شوند.

اجرای پیشنهادی از داخل `member2_kmer_analysis/`:

<div dir="ltr">

```bash
python scripts/kmer_analysis.py
python scripts/entropy.py
python scripts/similarity.py
python scripts/clustering.py
python scripts/tata_box.py
python scripts/stats_tests.py
python scripts/visualization.py
# شکل‌های آمار (در صورت نیاز جداگانه):
python scripts/stats_visualization.py
```

</div>

---

## ۶. معرفی فایل‌های ژنومی

| دیتاست | اسمبلی انتخاب‌شده | گونه (از data_summary.tsv) |
|---|---|---|
| `ncbi_dataset1` | `GCF_000002985.6_WBcel235_genomic.fna` | *Caenorhabditis elegans* |
| `ncbi_dataset2` | `GCF_000182925.2_NC12_genomic.fna` | *Neurospora crassa* |
| `ncbi_dataset3` | `GCF_001761485.1_ASM176148v1_genomic.fna` | *Yarrowia lipolytica* |
| `ncbi_dataset4` | `GCF_000004695.1_dicty_2.7_genomic.fna` | *Dictyostelium discoideum* |
| `ncbi_dataset5` | `GCA_000710315.1_Eremothecium_coryli_genomic.fna` | *Eremothecium coryli* |

همه کروموزوم‌ها/اسکفولدهای هر اسمبلی برای بردار k-mer گونه‌ای تجمیع شدند. برای آمار گروهی، همان رکوردها جداگانه نیز تحلیل شدند.

### ۶٫۱ ارجاع به جدول جمع‌آوری داده (بخش اول — فقط‌خواندنی)

پاسخ کامل پنج سؤال بخش اول پروژه (پایگاه دانلود، حجم فایل، تعداد کروموزوم و بازهای هر کروموزوم، نسخه اسمبلی، اندازه کل ژنوم) در خروجی دانشجوی ۱ تولید و مستند شده است و در این ماژول **تکرار نمی‌شود**. منبع مرجع فقط‌خواندنی:

- `member1_gc_analysis/results/tables/student1_gc_genome_summary.csv` (خلاصه سطح‌ژنوم: حجم فایل، accession/اسمبلی، تعداد کروموزوم/اسکفولد، اندازه کل، GC%)
- `member1_gc_analysis/results/tables/student1_gc_sequence_summary.csv` (طول باز به ازای هر رکورد/کروموزوم)

طبق گزارش دانشجوی ۱، برچسب پایگاه دانلود در محتوای فایل‌های ارائه‌شده صریحاً موجود نبود (`Not available in the provided file`)؛ سایر فیلدهای فوق از همان CSVها قابل خواندن‌اند.

---

## ۷. روش محاسبه K-mer

برای مقایسه اصلی و ماتریس شباهت/فاصله/PCA، k=۴ استفاده شد (بردار ۴⁴ = ۲۵۶ بُعد). علاوه بر آن، پرتکرارترین k-mer برای k=۳، k=۴ و k=۵ جداگانه محاسبه و در `results/kmer_top_by_k.csv` ذخیره شد. K-merهای گونه‌اختصاصی با تفاضل مجموعهٔ k-merهای پنج گونه (`results/species_unique_kmers.csv`) به‌دست آمد.

<div dir="ltr">

```text
f(kmer) = count(kmer) / Σ count(all kmers)
```

</div>

---

## ۸. روش Shannon Entropy

<div dir="ltr">

```text
H = −Σ pᵢ log₂(pᵢ)
```

</div>

که pᵢ فراوانی نسبی هر k-mer است. حداکثر نظری برای ۲۵۶ دسته برابر log₂(۲۵۶)=۸ بیت است.

---

## ۹. فاصله اقلیدسی و Cosine Similarity

اگر دو پروفایل نرمال‌شده u و v باشند:

- فاصله اقلیدسی: ‖u − v‖₂ (صفر یعنی پروفایل یکسان)
- cosine similarity: (u·v) / (‖u‖ ‖v‖) (یک یعنی جهت یکسان)

---

## ۱۰. PCA و خوشه‌بندی

ماتریس k-mer قبل از PCA با StandardScaler استانداردسازی شد. دو مؤلفه اصلی استخراج شدند. خوشه‌بندی سلسله‌مراتبی با پیوند Ward روی ماتریس فاصله اقلیدسی اجرا و دندروگرام رسم شد.

---

## ۱۱. جستجوی TATA Box

در ۵۰۰ bp اول توالی تجمیع‌شده هر گونه، موتیف متعارف `TATAAA` و واریانت‌های ساده‌تر جستجو شد. این اسکن توالی‌محور است و جایگزین annotation پروموترهای واقعی نیست.

---

## ۱۲. روش تحلیل آماری (بخش پنجم)

دو خانواده آزمون پیاده شد (`scripts/stats_tests.py`):

1. **همبستگی سطح‌گونه (n=۵):** Pearson، Spearman و permutation test (۹۹۹ resampling) برای  
   Genome Size ↔ Entropy / K-mer Diversity،  
   GC% ↔ Genome Size / Entropy / K-mer Diversity، و  
   TATA-box count ↔ GC%  
   (GC% فقط‌خواندنی از `member1_gc_analysis/results/tables/student1_gc_genome_summary.csv`).

2. **آزمون‌های گروهی روی رکوردها:** یک مقدار entropy و diversity برای هر کروموزوم/اسکفولد؛ سپس ANOVA یک‌طرفه، Kruskal–Wallis، permutation test روی آمارهٔ ANOVA، و Mann–Whitney زوجی با تصحیح Bonferroni.

**هشدار توان آماری:** با n=۵ در سطح گونه، همبستگی‌ها اکتشافی‌اند. در سطح رکورد، مشاهدات داخل یک ژنوم تکرار زیستی مستقل نیستند (nested / شبه-تکرار).

---

## ۱۳. نتایج خلاصه K-mer

| گونه | طول توالی معتبر (bp) | تعداد k-mer | انواع یکتا | پرتکرارترین (k=۴) | entropy (بیت) | kmer_diversity |
|---|---:|---:|---:|---|---:|---:|
| *C. elegans* | ۱۰۰٬۲۸۶٬۴۰۱ | ۱۰۰٬۲۸۶٬۳۹۸ | ۲۵۶ | AAAA | ۷٫۶۳۶ | ۲٫۵۵×۱۰⁻⁶ |
| *N. crassa* | ۴۱٬۰۶۱٬۶۰۳ | ۴۱٬۰۶۱٬۶۰۰ | ۲۵۶ | TTTT | ۷٫۹۶۵ | ۶٫۲۳×۱۰⁻⁶ |
| *Y. lipolytica* | ۲۰٬۵۰۰٬۰۶۶ | ۲۰٬۵۰۰٬۰۶۳ | ۲۵۶ | TTTT | ۷٫۹۴۶ | ۱٫۲۵×۱۰⁻⁵ |
| *D. discoideum* | ۳۴٬۱۸۱٬۸۳۱ | ۳۴٬۱۸۱٬۸۲۸ | ۲۵۶ | AAAA | ۶٫۹۱۷ | ۷٫۴۹×۱۰⁻⁶ |
| *E. coryli* | ۹٬۰۲۹٬۴۶۹ | ۹٬۰۲۹٬۴۶۶ | ۲۵۶ | AAAA | ۷٫۸۸۵ | ۲٫۸۴×۱۰⁻⁵ |

**پرتکرارترین k-mer برای k=۳، ۴ و ۵** (`kmer_top_by_k.csv`):

| گونه | k=۳ | k=۴ | k=۵ |
|---|---|---|---|
| *C. elegans* | TTT | AAAA | AAAAA |
| *N. crassa* | TTT | TTTT | TTTTT |
| *Y. lipolytica* | CAA | TTTT | AAAAA |
| *D. discoideum* | AAA | AAAA | AAAAA |
| *E. coryli* | AAA | AAAA | AAAAA |

**تعداد K-mer یکتا (k=۴):** در هر پنج گونه همهٔ ۲۵۶ نوع ۴-mer مشاهده شد.

**K-merهای گونه‌اختصاصی:** با تفاضل مجموعه‌ها در k=۳، ۴ و ۵، هیچ k-merی که فقط در یک گونه باشد یافت نشد (`species_unique_kmers.csv` خالی است). در این طول‌های کوتاه و ژنوم‌های بزرگ، الفبای ممکن تقریباً کامل در همه گونه‌ها دیده می‌شود؛ تمایز عمدتاً در فراوانی نسبی است، نه در حضور/غیاب.

**بیشترین تنوع K-mer (رتبه‌بندی بر اساس `kmer_diversity = unique/total`):** *E. coryli* (۲٫۸۴×۱۰⁻⁵) > *Y. lipolytica* > *D. discoideum* > *N. crassa* > *C. elegans*. توجه: این معیار با Shannon entropy یکی نیست؛ بالاترین entropy متعلق به *N. crassa* است.

**آیا K-merهای مشابه بیانگر نزدیکی تکاملی‌اند؟** در این نمونه، پروفایل‌های بسیار مشابه *N. crassa* و *Y. lipolytica* با ادغام زودهنگام آن‌ها در دندروگرام Ward و نزدیکی در فضای PCA هم‌خوان است؛ بااین‌حال درخت ترکیبی جایگزین فیلوژنی رسمی نیست و شباهت k-mer را باید به‌عنوان نزدیکی ترکیبی تفسیر کرد، نه اثبات قطعی خویشاوندی.

---

## ۱۴. نتایج Shannon Entropy

بیشترین entropy: *N. crassa* (۷٫۹۶۵) و سپس *Y. lipolytica* (۷٫۹۴۶). کمترین: *D. discoideum* (۶٫۹۱۷)، هم‌راستا با ترکیب AT-غنی و تسلط AAAA.

---

## ۱۵. ماتریس فاصله و شباهت

نزدیک‌ترین جفت: *N. crassa*–*Y. lipolytica* (فاصله اقلیدسی ۰٫۰۱۵۴؛ cosine ۰٫۹۷۱).  
دورترین جفت: *Y. lipolytica*–*D. discoideum* (فاصله اقلیدسی ۰٫۱۰۲۳؛ cosine ۰٫۵۴۲).  
*C. elegans* با *E. coryli* فاصله متوسط ≈ ۰٫۰۴۰ داشت.

**آیا گونه‌های جانوری در یک خوشه قرار می‌گیرند؟** با دادهٔ فعلی قابل پاسخ نیست: فقط یک گونهٔ جانوری (*C. elegans*) در مجموعه وجود دارد و نمی‌توان خوشه‌بندی «جانوران» را آزمود.  
**آیا گیاه از جانوران جدا می‌شود؟** با دادهٔ فعلی قابل پاسخ نیست: هیچ گیاهی در پنج اسمبلی موجود نیست (مجموعه شامل یک نماتود، سه قارچ و یک قارچ ژله‌ای/آمیبوزوآ است).  
**آیا نتایج با طبقه‌بندی زیستی مطابقت دارند؟** تا حدی بله در سطح توصیفی: دو قارچ *N. crassa* و *Y. lipolytica* نزدیک‌ترین جفت‌اند و *D. discoideum* (Amoebozoa، AT-rich) از قارچ‌ها جدا می‌شود؛ بااین‌حال این یک درخت فیلوژنتیک رسمی نیست و نباید فراتر از الگوی ترکیبی k-mer تعمیم داده شود.

---

## ۱۶. نتایج PCA و خوشه‌بندی

PC1 ≈ ۷۰٫۳۶٪ و PC2 ≈ ۱۲٫۴۴٪ واریانس (مجموع ≈ ۸۲٫۸٪). در فضای PCA، دو قارچ نزدیک هم و *D. discoideum* جدا بودند. دندروگرام Ward ابتدا دو قارچ را ادغام کرد، سپس *E. coryli*، بعد *C. elegans* و در نهایت *D. discoideum*.

---

## ۱۷. نتایج تحلیل آماری

### ۱۷٫۱ همبستگی سطح‌گونه (n=۵)

| متغیر X | متغیر Y | Pearson r | p | Spearman ρ | p | permutation p | تفسیر کوتاه |
|---|---|---:|---:|---:|---:|---:|---|
| Genome size | Entropy | −۰٫۱۵۱ | ۰٫۸۰۹ | −۰٫۱۰۰ | ۰٫۸۷۳ | ۰٫۶۳۶ | رابطه ضعیف؛ شواهد کافی نیست |
| Genome size | Diversity | −۰٫۷۴۳ | ۰٫۱۵۰ | −۱٫۰۰۰ | ≈۰ | ۰٫۰۱۸ | اکتشافی با n کوچک؛ با احتیاط |
| GC% | Genome size | −۰٫۲۴۹ | ۰٫۶۸۷ | −۰٫۳۰۰ | ۰٫۶۲۴ | ۰٫۶۰۰ | رابطه ضعیف و غیرمعنادار در این نمونه |
| GC% | Entropy | ۰٫۹۶۶ | ۰٫۰۰۷ | ۰٫۹۰۰ | ۰٫۰۳۷ | ۰٫۰۴۴ | همبستگی مثبت قوی در این نمونه پنج‌تایی |
| GC% | Diversity | ۰٫۲۳۳ | ۰٫۷۰۶ | ۰٫۳۰۰ | ۰٫۶۲۴ | ۰٫۷۸۸ | شواهد کافی نیست |
| TATA count | GC% | ۰٫۱۱۲ | ۰٫۸۵۸ | −۰٫۱۱۲ | ۰٫۸۵۸ | ۰٫۹۸۸ | شواهدی برای ارتباط TATA–GC در این نمونه نیست |

**پاسخ‌های مستقیم بخش پنجم:**

1. **آیا تفاوت GC بین گونه‌ها معنی‌دار است؟** در این ماژول آزمون گروهی روی GC اجرا نشد؛ طبق گزارش دانشجوی ۱ روی رکوردهای کروموزومی اصلی، ANOVA با F≈۶۱۶۴٫۲۴ و p≈۴٫۴۹×۱۰⁻³¹ تفاوت GC را معنادار نشان داد (با احتیاط شبه-تکرار).
2. **آیا تفاوت Shannon Entropy معنی‌دار است؟** بله در سطح رکورد: ANOVA F≈۱۹٫۸۱، p≈۱٫۰۹×۱۰⁻¹¹؛ Kruskal–Wallis H≈۴۲٫۷۶، p≈۱٫۱۶×۱۰⁻⁸؛ permutation p≈۰٫۰۰۱.
3. **آیا Genome Size با K-mer Diversity ارتباط دارد؟** در سطح گونه Pearson r≈−۰٫۷۴۳ (p≈۰٫۱۵۰) غیرمعنادار است؛ Spearman ρ≈−۱ و permutation p≈۰٫۰۱۸ اکتشافی‌اند و با n=۵ نباید قطعی ادعا شوند.
4. **آیا GC با پیچیدگی ژنوم مرتبط است؟** اگر پیچیدگی با Shannon entropy توزیع k-mer سنجیده شود، بله در این نمونه: Pearson r≈۰٫۹۶۶، p≈۰٫۰۰۷ (permutation p≈۰٫۰۴۴).

### ۱۷٫۲ آزمون‌های گروهی (روی رکوردها)

| متریک | ANOVA F | ANOVA p | Kruskal H | Kruskal p |
|---|---:|---:|---:|---:|
| Shannon entropy | ۱۹٫۸۱ | ≈۱٫۱×۱۰⁻¹¹ | ۴۲٫۷۶ | ≈۱٫۲×۱۰⁻⁸ |
| K-mer diversity | ۱۳٫۳۸ | ≈۱٫۴×۱۰⁻⁸ | ۳۹٫۶۷ | ≈۵٫۱×۱۰⁻⁸ |

Mann–Whitney زوجی با Bonferroni نشان داد بسیاری از جفت‌گونه‌ها از نظر entropy یا diversity تفاوت معنادار دارند؛ برای مثال *C. elegans* در برابر *D. discoideum* و *E. coryli* در entropy پس از تصحیح معنادار بود، درحالی‌که *N. crassa* در برابر *Y. lipolytica* در entropy پس از Bonferroni معنادار نبود (سازگار با شباهت ترکیبی بالای آن‌ها).

---

## ۱۸. نتایج TATA Box

| گونه | TATAAA در ۵۰۰ bp اول | تعداد | موقعیت(ها) |
|---|---|---:|---|
| *C. elegans* | خیر | ۰ | — |
| *N. crassa* | بله | ۵ | ۳۱۸، ۳۴۱، ۳۵۰، ۳۶۴، ۳۷۷ |
| *Y. lipolytica* | خیر | ۰ | — |
| *D. discoideum* | بله | ۲ | ۶۲، ۲۴۲ |
| *E. coryli* | خیر (فقط واریانت) | ۰ متعارف | — |

نبود TATAAA در پنجره کوتاه ≠ نبود پروموتر در کل ژنوم.

**پاسخ‌های مستقیم بخش ششم:**

1. **چه درصدی از ژن‌ها دارای TATA Box هستند؟** با داده/روش فعلی قابل پاسخ نیست: اسکن روی ۵۰۰ bp اول توالی ژنومی تجمیع‌شده است، نه روی نواحی پروموتر annotated ژن‌ها (GFF/GTF در دسترس نبود)؛ بنابراین نمی‌توان «درصد ژن‌ها» را گزارش کرد.
2. **آیا فراوانی TATA بین گونه‌ها متفاوت است؟** در همین پنجرهٔ محدود، بله توصیفی: *N. crassa* ۵ محل و *D. discoideum* ۲ محل داشتند و سه گونه دیگر صفر؛ آزمون آماری گروهی روی این شمارش‌های سطح‌گونه با n=۵ توان کافی ندارد.
3. **آیا ژن‌های دارای TATA بیان بیشتری دارند؟** با دادهٔ فعلی قابل پاسخ نیست: هیچ دادهٔ بیان ژن (RNA-seq/microarray) در پروژه جمع‌آوری نشده است.
4. **آیا گونه‌های نزدیک تکاملی الگوهای مشابهی دارند؟** در این نمونه الگو یکدست نیست: دو قارچ نزدیک (*N. crassa* و *Y. lipolytica*) یکی ۵ و دیگری ۰ محل TATAAA در پنجره ۵۰۰ bp داشتند؛ بنابراین شباهت تکاملی الزاماً الگوی TATA یکسان در این اسکن کوتاه ایجاد نکرد.
5. **آیا توزیع TATA با GC مرتبط است؟** در سطح گونه همبستگی ضعیف و غیرمعنادار بود (Pearson r≈۰٫۱۱۲، p≈۰٫۸۵۸؛ Spearman ρ≈−۰٫۱۱۲، p≈۰٫۸۵۸؛ permutation p≈۰٫۹۸۸).

---

## ۱۹. تفسیر نمودارها

**شکل‌های اصلی K-mer**

1. `kmer_frequency_heatmap.png` — فراوانی نسبی k-merها در میان گونه‌ها (نه فاصله)  
2. `heatmap.png` — فاصله اقلیدسی زوجی  
3. `similarity_matrix.png` — cosine similarity (بخش چهارم)  
4. `pca.png` — پراکنش در فضای PC1–PC2  
5. `dendrogram.png` — درخت خوشه‌بندی Ward  
6. `entropy_comparison.png` — مقایسه Shannon entropy گونه‌ها  
7. `correlation_matrix.png` — ماتریس همبستگی بخش پنجم (GC، اندازه، entropy، diversity)  

**شکل‌های آماری تکمیلی**

8. `correlation_scatter.png` — پراکنش همبستگی‌های سطح‌گونه  
9. `entropy_by_organism_boxplot.png` — پراکندگی entropy بین رکوردهای هر گونه  
10. `diversity_by_organism_boxplot.png` — پراکندگی k-mer diversity بین رکوردها  
11. `anova_kruskal_summary.png` — خلاصه ANOVA / Kruskal–Wallis  
12. `mannwhitney_heatmap.png` — معناداری زوجی Mann–Whitney (پس از Bonferroni)

---

## ۲۰. بحث

پروفایل k-mer ساختار روابط پنج ژنوم را بدون هم‌ترازسازی کامل نشان داد. شباهت بالای دو قارچ و جدایی *D. discoideum* با ترکیب AT-غنی و entropy پایین‌تر سازگار است. همبستگی مثبت GC با entropy در این نمونه کوچک با الگوی «ژنوم AT-rich با توزیع نایکنواخت‌تر k-mer» هم‌خوان است، ولی نباید به قانون عمومی تعمیم داده شود. آزمون‌های گروهی تفاوت بین گونه‌ها را تقویت می‌کنند، اما به‌خاطر nested بودن رکوردها باید محتاطانه تفسیر شوند. تحلیل TATA box محدود به پنجره ۵۰۰ bp و ترتیب رکوردهای تجمیع‌شده است.

---

## ۲۱. محدودیت‌ها

1. k=۴ تفکیک ریزمقیاس توالی‌های بسیار نزدیک را محدود می‌کند.  
2. تجمیع contigها اثر ترکیب محلی (کروموزوم در برابر اندامک) را رقیق می‌کند.  
3. با n=۵، خوشه‌بندی، PCA و همبستگی سطح‌گونه توصیفی‌اند.  
4. رکوردهای یک ژنوم تکرار زیستی مستقل نیستند.  
5. TATA box فقط در ۵۰۰ bp اول توالی تجمیع‌شده جستجو شد.  
6. نتایج به کیفیت اسمبلی و نسخه RefSeq/GenBank حساس‌اند.  
7. تحلیل GC متعلق به دانشجوی ۱ است؛ اینجا فقط برای همبستگی اختیاری خوانده می‌شود.

---

## ۲۲. پاسخ مستقیم به سؤالات بخش K-mer و آمار

- **فراوان‌ترین 3-mer / 4-mer / 5-mer؟** جدول §۱۳ (`kmer_top_by_k.csv`).  
- **تعداد K-mer یکتا؟** در k=۴ همه گونه‌ها ۲۵۶ نوع.  
- **کدام K-mer فقط در یک گونه؟** در k=۳–۵ هیچ‌کدام (مجموعه خالی).  
- **بیشترین تنوع K-mer؟** *E. coryli* بر اساس `kmer_diversity`.  
- **نزدیک‌ترین / دورترین گونه‌ها؟** نزدیک: *N. crassa*–*Y. lipolytica*؛ دور: *Y. lipolytica*–*D. discoideum*.  
- **خوشه جانوران / جدایی گیاه؟** با دادهٔ فعلی غیرقابل‌پاسخ (فقط یک جانور؛ بدون گیاه).  
- **آیا تفاوت entropy بین گونه‌ها معنادار است؟** بله در سطح رکورد (ANOVA/KW؛ §۱۷٫۲).  
- **آیا Genome Size با Diversity ارتباط دارد؟** Pearson غیرمعنادار؛ نتایج اکتشافی با n=۵.  
- **آیا GC با پیچیدگی مرتبط است؟** با entropy به‌عنوان پروکسی، بله در این نمونه (r≈۰٫۹۶۶، p≈۰٫۰۰۷).  
- **TATA: درصد ژن‌ها / بیان / گیاه–جانور؟** سه مورد غیرقابل‌پاسخ به‌خاطر نبود annotation ژن، داده بیان، و گیاه در مجموعه (§۱۸).  
- **TATA با GC؟** همبستگی غیرمعنادار (r≈۰٫۱۱۲، p≈۰٫۸۵۸).

---

## ۲۳. نتیجه‌گیری

پایپ‌لاین دانشجوی ۲ روی هر پنج دیتاست NCBI اجرا شد و خروجی‌های جدولی، آماری و تصویری قابل بازتولید تولید کرد. پروفایل ۴-mer، entropy، ماتریس‌های شباهت/فاصله، PCA، دندروگرام، اسکن TATA box و آزمون‌های آماری تصویر سازگاری از روابط ترکیبی این پنج ژنوم ارائه دادند. مناسب‌ترین تفسیر علمی، توصیفی و محدود به همین پنج اسمبلی است.

---

## ۲۴. فهرست فایل‌های خروجی

**اسکریپت‌ها**

- `scripts/kmer_analysis.py`
- `scripts/entropy.py`
- `scripts/similarity.py`
- `scripts/clustering.py`
- `scripts/tata_box.py`
- `scripts/stats_tests.py`
- `scripts/stats_visualization.py`
- `scripts/visualization.py`

**نتایج اصلی و آمار**

- `results/kmer_summary.csv`
- `results/kmer_top_by_k.csv`
- `results/species_unique_kmers.csv`
- `results/entropy.csv`
- `results/distance_matrix.csv`
- `results/cosine_similarity.csv`
- `results/tata_box.csv`
- `results/pca_coordinates.csv`
- `results/pca_variance.csv`
- `results/linkage_matrix.csv`
- `results/species_level_metrics.csv`
- `results/per_record_metrics.csv`
- `results/correlation_results.csv`
- `results/anova_kruskal_results.csv`
- `results/mannwhitney_pairwise.csv`

**شکل‌ها**

- `results/figures/kmer_frequency_heatmap.png`
- `results/figures/heatmap.png`
- `results/figures/similarity_matrix.png`
- `results/figures/pca.png`
- `results/figures/dendrogram.png`
- `results/figures/entropy_comparison.png`
- `results/figures/correlation_matrix.png`
- `results/figures/correlation_scatter.png`
- `results/figures/entropy_by_organism_boxplot.png`
- `results/figures/diversity_by_organism_boxplot.png`
- `results/figures/anova_kruskal_summary.png`
- `results/figures/mannwhitney_heatmap.png`

**مستندات**

- `report/student2_kmer_report_fa.md`
- `report/student2_kmer_oral_presentation_guide_fa.md`
- `README.md`
- `requirements.txt`

</div>
