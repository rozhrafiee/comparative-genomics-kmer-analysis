<div dir="rtl" lang="fa">

# ۱. عنوان پروژه

**ژنومیک مقایسه‌ای: تحلیل K-mer و محتوای GC در پنج ارگانیسم — بخش دانشجوی ۲**

این گزارش فقط بخش تحلیل K-mer، Shannon Entropy، ماتریس شباهت و فاصله، PCA، خوشه‌بندی سلسله‌مراتبی، دندروگرام، جستجوی TATA box، تحلیل آماری (همبستگی، ANOVA، Kruskal–Wallis، Mann–Whitney) و نمودارهای مربوط به همین تحلیل‌ها را پوشش می‌دهد.

---

## ۲. چکیده

پنج اسمبلی ژنومی NCBI (`ncbi_dataset1` تا `ncbi_dataset5`) با پایپ‌لاین قابل بازتولید در پوشه `member2_kmer_analysis` تحلیل شدند. هویت گونه‌ها از `data_summary.tsv` استخراج شد: *Caenorhabditis elegans*، *Neurospora crassa*، *Yarrowia lipolytica*، *Dictyostelium discoideum* و *Eremothecium coryli*. فهرست گونه‌های Brief فقط نمونه بود؛ تحلیل روی همین پنج اسمبلی واقعی انجام شد.

برای هر گونه، همه کروموزوم‌ها/اسکفولدها در یک پروفایل k-mer تجمیع شدند (k=۴). هر پنج ژنوم همه ۲۵۶ نوع ۴-mer ممکن را داشتند. بیشترین Shannon entropy مربوط به *N. crassa* (۷٫۹۶۵ بیت) و کمترین مربوط به *D. discoideum* (۶٫۹۱۷ بیت) بود. نزدیک‌ترین جفت از نظر فاصله اقلیدسی *N. crassa* و *Y. lipolytica* (۰٫۰۱۵۴) و cosine similarity آن‌ها ۰٫۹۷۱ بود. در PCA، PC1 حدود ۷۰٫۴٪ واریانس را توضیح داد.

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

برای آزمون‌های گروهی (ANOVA / Kruskal–Wallis / Mann–Whitney)، تابع `load_sequences_by_record` هر رکورد FASTA را جدا نگه می‌دارد تا برای هر گونه چند نمونه (کروموزوم/اسکفولد) وجود داشته باشد. خروجی‌های آمار در `results/statistics/` و شکل‌های آمار در `results/figures/` ذخیره می‌شوند.

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

---

## ۷. روش محاسبه K-mer

با k=۴، پنجره لغزنده روی توالی حرکت می‌کند و هر زیررشته ۴ نوکلئوتیدی شمارش می‌شود. بردار نهایی دارای ۴⁴ = ۲۵۶ بُعد است.

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

1. **همبستگی سطح‌گونه (n=۵):** Pearson و Spearman برای  
   Genome Size ↔ Entropy / K-mer Diversity و  
   GC% ↔ Entropy / K-mer Diversity  
   (GC% فقط‌خواندنی از `member1_gc_analysis/results/tables/student1_gc_genome_summary.csv`).

2. **آزمون‌های گروهی روی رکوردها:** یک مقدار entropy و diversity برای هر کروموزوم/اسکفولد؛ سپس ANOVA یک‌طرفه، Kruskal–Wallis، و Mann–Whitney زوجی با تصحیح Bonferroni.

**هشدار توان آماری:** با n=۵ در سطح گونه، همبستگی‌ها اکتشافی‌اند. در سطح رکورد، مشاهدات داخل یک ژنوم تکرار زیستی مستقل نیستند (nested / شبه-تکرار).

---

## ۱۳. نتایج خلاصه K-mer

| گونه | طول توالی معتبر (bp) | تعداد k-mer | انواع یکتا | پرتکرارترین | entropy (بیت) |
|---|---:|---:|---:|---|---:|
| *C. elegans* | ۱۰۰٬۲۸۶٬۴۰۱ | ۱۰۰٬۲۸۶٬۳۹۸ | ۲۵۶ | AAAA | ۷٫۶۳۶ |
| *N. crassa* | ۴۱٬۰۶۱٬۶۰۳ | ۴۱٬۰۶۱٬۶۰۰ | ۲۵۶ | TTTT | ۷٫۹۶۵ |
| *Y. lipolytica* | ۲۰٬۵۰۰٬۰۶۶ | ۲۰٬۵۰۰٬۰۶۳ | ۲۵۶ | TTTT | ۷٫۹۴۶ |
| *D. discoideum* | ۳۴٬۱۸۱٬۸۳۱ | ۳۴٬۱۸۱٬۸۲۸ | ۲۵۶ | AAAA | ۶٫۹۱۷ |
| *E. coryli* | ۹٬۰۲۹٬۴۶۹ | ۹٬۰۲۹٬۴۶۶ | ۲۵۶ | AAAA | ۷٫۸۸۵ |

تفاوت‌ها عمدتاً در فراوانی نسبی است، نه در حضور/غیاب کامل ۴-merها.

---

## ۱۴. نتایج Shannon Entropy

بیشترین entropy: *N. crassa* (۷٫۹۶۵) و سپس *Y. lipolytica* (۷٫۹۴۶). کمترین: *D. discoideum* (۶٫۹۱۷)، هم‌راستا با ترکیب AT-غنی و تسلط AAAA.

---

## ۱۵. ماتریس فاصله و شباهت

نزدیک‌ترین جفت: *N. crassa*–*Y. lipolytica* (فاصله ۰٫۰۱۵۴؛ cosine ۰٫۹۷۱).  
دورترین الگوها: *D. discoideum* نسبت به قارچ‌ها (مثلاً فاصله با *Y. lipolytica* ≈ ۰٫۱۰۲؛ cosine ≈ ۰٫۵۴۲).  
*C. elegans* با *E. coryli* فاصله متوسط ≈ ۰٫۰۴۰ داشت.

---

## ۱۶. نتایج PCA و خوشه‌بندی

PC1 ≈ ۷۰٫۳۶٪ و PC2 ≈ ۱۲٫۴۴٪ واریانس (مجموع ≈ ۸۲٫۸٪). در فضای PCA، دو قارچ نزدیک هم و *D. discoideum* جدا بودند. دندروگرام Ward ابتدا دو قارچ را ادغام کرد، سپس *E. coryli*، بعد *C. elegans* و در نهایت *D. discoideum*.

---

## ۱۷. نتایج تحلیل آماری

### ۱۷٫۱ همبستگی سطح‌گونه (n=۵)

| متغیر X | متغیر Y | Pearson r | p | Spearman ρ | p | تفسیر کوتاه |
|---|---|---:|---:|---:|---:|---|
| Genome size | Entropy | −۰٫۱۵۱ | ۰٫۸۰۹ | −۰٫۱۰۰ | ۰٫۸۷۳ | رابطه ضعیف؛ شواهد کافی نیست |
| Genome size | Diversity | −۰٫۷۴۳ | ۰٫۱۵۰ | −۱٫۰۰۰ | ≈۰ | اکتشافی با n کوچک؛ با احتیاط |
| GC% | Entropy | ۰٫۹۶۶ | ۰٫۰۰۷ | ۰٫۹۰۰ | ۰٫۰۳۷ | همبستگی مثبت قوی در این نمونه پنج‌تایی |
| GC% | Diversity | ۰٫۲۳۳ | ۰٫۷۰۶ | ۰٫۳۰۰ | ۰٫۶۲۴ | شواهد کافی نیست |

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

---

## ۱۹. تفسیر نمودارها

**شکل‌های اصلی K-mer**

1. `heatmap.png` — فاصله اقلیدسی زوجی  
2. `similarity_matrix.png` — cosine similarity  
3. `pca.png` — پراکنش در فضای PC1–PC2  
4. `dendrogram.png` — درخت خوشه‌بندی Ward  
5. `entropy_comparison.png` — مقایسه Shannon entropy گونه‌ها  

**شکل‌های آماری جدید**

6. `correlation_scatter.png` — پراکنش همبستگی‌های سطح‌گونه  
7. `entropy_by_organism_boxplot.png` — پراکندگی entropy بین رکوردهای هر گونه  
8. `diversity_by_organism_boxplot.png` — پراکندگی k-mer diversity بین رکوردها  
9. `anova_kruskal_summary.png` — خلاصه ANOVA / Kruskal–Wallis  
10. `mannwhitney_heatmap.png` — معناداری زوجی Mann–Whitney (پس از Bonferroni)

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

- **K-mer چیست؟** زیررشته طول ثابت k از DNA.  
- **چرا نرمال‌سازی؟** برای مقایسه ژنوم‌های با طول متفاوت.  
- **Entropy چه می‌گوید؟** یکنواختی توزیع k-mer؛ نزدیک به ۸ یعنی پخش‌شدگی زیاد.  
- **کدام گونه‌ها شبیه‌ترند؟** *N. crassa* و *Y. lipolytica*.  
- **کدام جداست؟** *D. discoideum*.  
- **آیا تفاوت entropy بین گونه‌ها معنادار است؟** در سطح رکورد بله (ANOVA/KW)، با احتیاط شبه-تکرار.  
- **آیا Genome Size با Diversity ارتباط دارد؟** در n=۵ اکتشافی است و نباید قطعی ادعا شود.  
- **آیا GC با پیچیدگی مرتبط است؟** در این نمونه، GC با Shannon entropy همبستگی مثبت قوی نشان داد.  
- **TATA box کجا دیده شد؟** در پنجره ۵۰۰ bp برای *N. crassa* و *D. discoideum*.

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

- `results/kmer_summary.csv` (یا نسخه هم‌تراز در `results/statistics/`)
- `results/entropy.csv`
- `results/distance_matrix.csv`
- `results/cosine_similarity.csv`
- `results/tata_box.csv`
- `results/pca_coordinates.csv`
- `results/pca_variance.csv`
- `results/linkage_matrix.csv`
- `results/statistics/species_level_metrics.csv`
- `results/statistics/per_record_metrics.csv`
- `results/statistics/correlation_results.csv`
- `results/statistics/anova_kruskal_results.csv`
- `results/statistics/mannwhitney_pairwise.csv`

**شکل‌ها**

- `results/figures/heatmap.png`
- `results/figures/similarity_matrix.png`
- `results/figures/pca.png`
- `results/figures/dendrogram.png`
- `results/figures/entropy_comparison.png`
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
