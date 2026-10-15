"""Chinese revision and submission memo, including verified publisher links."""
from pathlib import Path
import json
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.platypus import SimpleDocTemplate,Paragraph,Spacer,PageBreak
from reportlab.lib.styles import ParagraphStyle
from reportlab.lib.pagesizes import A4
from reportlab.lib.units import cm
from reportlab.lib.enums import TA_LEFT

ROOT=Path(__file__).resolve().parent
SECTIONS=[
('修订结果与投稿建议','核对日期：2026年10月15日。作者、单位与通讯信息按要求保留占位；已注明无专项基金支持、无利益冲突。原 paper_v2 和原始结果保留，修订输出位于 paper_v3 与 results_v3。'),
('本轮完成的修改','重新绘制图1，区分真实模拟与示意网络；重排图4的图例、分支与势能曲线，将折叠细节移到独立补充图；重排图6，拆开颜色条、图例与数据区，移除遮挡曲线的文字。15个展示公式采用统一字号、分行和连续编号；正文参数下标与稳态上标同步校正。'),
('影响结论的实质性修正','图4的灰色缺口来自求根失败，不是无准稳态解的证据。新增 q=1 解析消元与括区间求解，准稳态残差最大约3.0×10⁻¹¹。定位两个折叠点0.928210与0.990198；迟滞区间宽0.061988。将双稳态轨迹延长到2000，低损伤平衡约0.0202，而120时刻约0.009尚未收敛。'),
('统计与引用修正','重算100个组合剂量点的精确匹配单药对照。以10⁻⁶容差排除求解器误差后，90%网格点为正Bliss超额，均值0.1529；基准清除率那一列为数值零。撤下原不成立的profile-likelihood置信区间，改报局部高斯不确定性并说明不能声称95%覆盖率。修复PRCC自由度。37篇实际引用的文献逐条核对DOI、作者与卷页，采用数字顺序Vancouver体例，并提供BibTeX与核对记录。'),
('当前科学定位','稿件应定位为机制假设与数值分析研究。没有患者数据、真实连接组或独立实验验证；二级成核项没有配套单体消耗，因此变量是现象学负荷而非严格质量守恒浓度。0.547是冻结单体池近似阈值，不是完整系统的精确阈值或普适治疗标准。200时刻的低损伤不等于已证明网络渐近停滞。已同步收紧摘要、结果与讨论中的过度外推。'),
('仍需增强的证据','若以较高质量SCI论文为目标，优先补充：一组公开或自有动力学数据的校准与外部验证；更符合质量守恒的模型对照；真实连接组或多随机图重复；与噪声模型一致的联合参数剖面与重复噪声覆盖率实验。这些证据不能用文字润色代替，本轮未虚构实验或临床结果。'),
('投稿顺序与费用选择','建议先考虑 Bulletin of Mathematical Biology；若进一步补强生物学证据，可考虑 Journal of Theoretical Biology；若补强解析证明，可尝试 Journal of Mathematical Biology。另可考虑 Mathematical Medicine and Biology，主题较匹配，但仍需核实是否适合当前无外部数据的证据强度。以下为期刊适配判断，不是录用承诺。'),
('1  Bulletin of Mathematical Biology','官方明确为SCIE收录的混合期刊，覆盖数学与生命科学交叉研究；选择订阅出版路径不收APC。本稿的稳态、分岔和干预机制较符合其定位。当前指南要求作者—年份引文与字母序参考文献，因此决定投此刊后需从当前数字顺序稿切换体例，并按要求加入连续行号。<br/><link href="https://link.springer.com/journal/11538">期刊与收录信息</link>；<link href="https://link.springer.com/journal/11538/how-to-publish-with-us">出版路径与费用</link>；<link href="https://link.springer.com/journal/11538/submission-guidelines">作者指南</link>。'),
('2  Journal of Theoretical Biology','范围包括脑与神经科学、网络与复杂系统及理论建模，强调新的生物学理解。现象学模型需更紧密地连接实验依据。本刊支持订阅出版；Elsevier官方政策明确混合期刊可选择免费的订阅路径。投稿前确认印刷彩色等可选项目，不选择付费OA。<br/><link href="https://shop.elsevier.com/journals/journal-of-theoretical-biology/0022-5193">官方范围与订阅信息</link>；<link href="https://www.elsevier.com/about/policies-and-standards/pricing">Elsevier费用政策</link>。'),
('3  Journal of Mathematical Biology','官方支持订阅出版，订阅路径不收APC。强调通过严格数学分析得到新的生物学认识；适合作为补强数学理论后的备选，当前稿主要依赖数值分析，不能仅因具有分岔图就视为充分匹配。<br/><link href="https://link.springer.com/journal/285/how-to-publish-with-us">费用政策</link>；<link href="https://link.springer.com/journal/285/aims-and-scope">范围与要求</link>。'),
('4  Mathematical Medicine and Biology','已发表AD聚集、扩散与疾病进展模型，主题匹配。官方提供标准许可和OA许可两条路径；本预算目标应选择标准许可，OA为另收费选项。作者指南采用Harvard作者—年份参考文献格式。<br/><link href="https://academic.oup.com/imammb/pages/General_Instructions">官方作者指南与许可选择</link>。'),
('交付与核验范围','主PDF已逐页渲染检查。Word副本保留可编辑正文与高分辨率公式图片；自动审批拒绝了Word排版核验所需的LibreOffice安装包解包命令，返回“blocked by policy”，因此未声称Word副本已完成视觉核验。数值回归检查覆盖快平衡、双稳态、传输守恒权重、PRCC自由度、文献与图片完整性及Bliss对照。'),
]

def main():
    out=ROOT/'paper_v3';out.mkdir(exist_ok=True)
    pdfmetrics.registerFont(TTFont('CN',r'C:/Windows/Fonts/simsun.ttc',subfontIndex=0))
    h=ParagraphStyle('h',fontName='CN',fontSize=12,leading=17,spaceBefore=12,spaceAfter=6,keepWithNext=True)
    p=ParagraphStyle('p',fontName='CN',fontSize=10.5,leading=17,spaceAfter=7,wordWrap='CJK')
    title=ParagraphStyle('title',parent=h,fontSize=17,leading=24,spaceAfter=12)
    flow=[]
    for i,(head,text) in enumerate(SECTIONS):
        text=text.replace('10⁻¹¹','10<super>-11</super>').replace('10⁻⁶','10<super>-6</super>')
        if i==6:flow.append(PageBreak())
        flow.append(Paragraph(head,title if i==0 else h));flow.append(Paragraph(text,p))
    def number(canvas,doc):
        canvas.setFont('CN',9);canvas.drawCentredString(A4[0]/2,1.2*cm,str(doc.page))
    SimpleDocTemplate(str(out/'Revision_and_Submission_Notes.pdf'),pagesize=A4,
        leftMargin=2.3*cm,rightMargin=2.3*cm,topMargin=2.0*cm,bottomMargin=2.0*cm).build(flow,onFirstPage=number,onLaterPages=number)
    import re
    (out/'REVISION_NOTES.txt').write_text('\n\n'.join(a+'\n'+re.sub('<[^>]+>','',b) for a,b in SECTIONS),encoding='utf-8')
    print('Revision and submission memo saved')

if __name__=='__main__':main()
