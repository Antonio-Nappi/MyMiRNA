import { Component, OnInit } from '@angular/core';
import { NgForm } from '@angular/forms';
import { RequestService } from './requests.service';
import { DomSanitizer, SafeResourceUrl } from '@angular/platform-browser';

@Component({
  selector: 'app-home',
  templateUrl: 'home.page.html',
  styleUrls: ['home.page.scss'],
})
export class HomePage implements OnInit {

  // Booleans that represents the current process. They have been used to show the relative ion-card
  public default: boolean;
  public isTrimming: boolean;
  public isShortStack: boolean;
  public isMiRNA: boolean;
  public isPiRNA: boolean;
  public isOthers: boolean;
  public isNovelMirna: boolean;
  public isNovelPirna: boolean;

  // HTML pages that represent the results of the relative analysis
  public multiqc: SafeResourceUrl;
  public shortStack: SafeResourceUrl;
  public mirnaAnalysis: SafeResourceUrl;
  public matureMirnaList: string[];
  public preMirnaList: string[];
  public mirnaInformation: SafeResourceUrl;

  constructor(private requestService: RequestService, private sanitizer: DomSanitizer) {
    this.default = false;
    this.isTrimming = false;
    this.isShortStack = false;
    this.isPiRNA = false;
    this.isOthers = false;
    this.isMiRNA = false;
    this.isNovelMirna = false;
    this.isNovelPirna = false;
  }

  ngOnInit() {
    document.getElementById('card2').setAttribute('disabled', 'true');
    document.getElementById('card3').setAttribute('disabled', 'true');
    document.getElementById('card4').setAttribute('disabled', 'true');
    document.getElementById('card5').setAttribute('disabled', 'true');
    document.getElementById('card6').setAttribute('disabled', 'true');
    document.getElementById('card7').setAttribute('disabled', 'true');
  }

  private buildBody(form: NgForm) {
    return JSON.stringify({
      multimap: form.value['mismatches'],
      cores: form.value['core'],
      p_value: form.value['p_value'],
      p_value_adjusted: form.value['p_value_adjusted'],
      log_2_fold: form.value['log']
    });
  }

  onSubmit(form: NgForm, name: string) {
    if (form.valid) {
      if (name === 'f1') {  // adapter trimming step
        const body = JSON.stringify({
          qual: form.value['qual'],
          adapter: form.value['adapter']
        });
        this.requestService.trimmingStep(body).subscribe(res => {
          this.multiqc = this.sanitizer.bypassSecurityTrustResourceUrl(res);
          document.getElementById('card2').setAttribute('disabled', 'false');
          this.isTrimming = true;
          this.default = true;
        });
      } else if (name === 'f2') { // shortstack step
        const body = JSON.stringify({
          multimap: form.value['tresh'],
          cores: form.value['core'],
          p_value: form.value['p_value'],
          p_value_adjusted: form.value['p_value_adjusted'],
          log_2_fold: form.value['log']
        });
        this.requestService.shortStack(body).subscribe(res => {
          this.shortStack = this.sanitizer.bypassSecurityTrustResourceUrl(res);
          document.getElementById('card3').setAttribute('disabled', 'false');
          this.isShortStack = true;
        });
      } else if (form.name === 'f3') {  // mirna analysis
        const body = this.buildBody(form);
        this.requestService.mirnaAnalysis(body).subscribe(res => {  // differential analysis
          this.mirnaAnalysis = this.sanitizer.bypassSecurityTrustResourceUrl(res);
          this.requestService.mirnaDetection().subscribe(mirnaList => { // mirna & pre-mirna information
            this.matureMirnaList = mirnaList.mature;
            this.preMirnaList = mirnaList.pre;
          });
          document.getElementById('card4').setAttribute('disabled', 'false');
          this.isMiRNA = true;
        });
      } else if (name === 'f4') { // pirna analysis
        const body = this.buildBody(form);
        // richiesta
        document.getElementById('card5').setAttribute('disabled', 'false');
      } else if (name === 'f5') {
        const body = this.buildBody(form);
        // richiesta
        document.getElementById('card6').setAttribute('disabled', 'false');
      } else if (name === 'f6') {
        const body = this.buildBody(form);
        // richiesta
        document.getElementById('card7').setAttribute('disabled', 'false');
      } else {
        const body = this.buildBody(form);
        // richiesta
      }
    }
  }

  showMature(mirna: string) {
    const path = 'http://localhost:8080/mirnas/' + mirna;
    this.requestService.mirnaInformation(path).subscribe(res => {
      this.mirnaInformation = this.sanitizer.bypassSecurityTrustResourceUrl(res);
      // ADD MODAL
    });
  }

}
